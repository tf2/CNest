/*
 * =====================================================================================
 *
 *       Filename:  cnv_baf_data.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  09/05/17 21:39:29
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
/* This program to extract copy number data from BAM / CRAM files
 * Author: Tomas William Fitzgerald
 * Email: tf2@sanger.ac.uk
 * To compile this program, you may:
 * Nb. Assumes you are within a samtools directory containing htslib-1.3
 *
 *   gcc -O3 -Wall -Isrc -I. -Ihtslib-1.3/ -Ihtslib-1.3/htslib -rdynamic -o cnv_baf_data -Lhtslib-1.3/ -L. cnv_baf_data.c -lhts -lbam -lpthread -lz -lm
 *
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include "htslib-1.3/htslib/sam.h" 
#include "samtools.h"

// Some fixed cut-offs used for the B allele (baf) calculation
const int 	pileup_depth_min_threshold	= 7;
const float min_het_baf                	= 0.06;

// Auxiliary data structure
typedef struct {     // file parsing information
    htsFile *fp;     // the file handle
    bam_hdr_t *hdr;  // the file header
    hts_itr_t *iter; // NULL if a region not specified
    hts_idx_t *idx; // the hts index
    int min_mapQ;
} fdata;

// Initialize auxiliary data structure
fdata *get_fdata (char *filename, int mapQ) {
	fdata *data;
	data = calloc(1, sizeof(fdata));
    data->fp = hts_open(filename, "r"); // open BAM
    if (hts_set_opt(data->fp, CRAM_OPT_REQUIRED_FIELDS,
                    SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR |
                    SAM_SEQ)) {
        fprintf(stderr, "Failed to set CRAM_OPT_REQUIRED_FIELDS value\n");
    }
    if (hts_set_opt(data->fp, CRAM_OPT_DECODE_MD, 0)) {
        fprintf(stderr, "Failed to set CRAM_OPT_DECODE_MD value\n");
    }
    data->min_mapQ = mapQ;                    // set the mapQ filter
    data->hdr = sam_hdr_read(data->fp);    // read the BAM header
    data->idx = sam_index_load(data->fp, filename);  // load the index
return data;
}

// This function reads a BAM alignment from one BAM file
static int read_bam(void *data, bam1_t *b) {
    fdata *aux = (fdata*)data; // data in fact is a pointer to an auxiliary structure
    int ret;
    while (1) {
        ret = aux->iter? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b);
        if ( ret<0 ) break;
        if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;
        if ( (int)b->core.qual < aux->min_mapQ ) continue;
        break;
    }
return ret;
}

// This function defines the genotype using a fixed baf cut-off
char* convert_baf_to_genotype (float baf) {
 	char *gt = "AB";
	if (baf < min_het_baf) gt="AA";
	if (baf > (1-min_het_baf)) gt="BB";
	if (baf == -1) gt="NA";
return gt;
}

// This function extracts the allele counts at a position and calculates the baf
float extract_allele_counts (fdata **data, char *reg, char *ref_allele) {
	int i,base_a =0,base_c=0,base_g=0,base_t=0,base_n=0;
    data[0]->iter = sam_itr_querys(data[0]->idx, data[0]->hdr, reg); // set the iterator
	int ret, n=1, tid, beg, end, pos, *n_plp;
	const bam_pileup1_t **plp; bam_mplp_t mplp;
    mplp = bam_mplp_init(n, read_bam, (void**)data); // initialization
    n_plp = calloc(n, sizeof(int)); // n_plp[i] is the number of covering reads from the i-th BAM
    plp = calloc(n, sizeof(bam_pileup1_t*)); 
    beg = data[0]->iter->beg; end = data[0]->iter->end;
    while ((ret=bam_mplp_auto(mplp, &tid, &pos, n_plp, plp)) > 0) { 
		if (pos < beg || pos >= end) continue;
		for(i=0;i<n_plp[0];i++) {
			const bam_pileup1_t *p = plp[0] + i;
			int8_t base = bam_seqi(bam_get_seq(p->b), p->qpos);
			switch ( base ) {
				case 1: base_a++; break;
				case 2: base_c++; break;
				case 4: base_g++; break;
				case 8: base_t++; break;
				case 15: base_n++; break;
			}
		}
		
    }    
    float baf = -1;
    int read_count = base_a+base_c+base_g+base_t+base_n;
    if(read_count>pileup_depth_min_threshold) {
    	switch ( ref_allele[0] ) {
			case 'A': baf = (float)(base_c+base_g+base_t+base_n)/(float)read_count; break;
			case 'C': baf = (float)(base_a+base_g+base_t+base_n)/(float)read_count; break;
			case 'G': baf = (float)(base_a+base_c+base_t+base_n)/(float)read_count; break;
			case 'T': baf = (float)(base_a+base_c+base_g+base_n)/(float)read_count; break;
			case 'N': baf = (float)(base_a+base_c+base_g+base_t)/(float)read_count; break;
		}	
	}
   free(n_plp); free(plp); 
   bam_mplp_destroy(mplp);
   sam_itr_destroy(data[0]->iter);
return baf;
}

// This function calculates the read count and depth across an interval
int read_cnv_data (fdata *data, char *reg, float *depth) {
	int i, result, read_count=0, coverage = 0;
	bam1_t* b = NULL; b = bam_init1();
	data->iter = sam_itr_querys(data->idx, data->hdr, reg); // set the iterator	
	while ((result = sam_itr_next(data->fp, data->iter, b)) >= 0) {
    	if(b->core.qual < data->min_mapQ || (b->core.flag & BAM_FUNMAP)
			|| !(b->core.flag & BAM_FPROPER_PAIR) || (b->core.flag & BAM_FMUNMAP)//Proper pair and mate unmapped
			|| (b->core.flag & BAM_FDUP)//1024 is PCR/optical duplicate
			|| (b->core.flag & BAM_FSECONDARY) || (b->core.flag & BAM_FQCFAIL)//Secondary alignment, quality fail
			|| (b->core.flag & BAM_FSUPPLEMENTARY) ) continue;
    	for(i=0;i<b->core.l_qseq;i++) {
    		if(b->core.pos >= data->iter->beg && b->core.pos+i <= data->iter->end) coverage += 1;
    	}
    	read_count++;
	}
	*depth = (float)coverage /(data->iter->end-data->iter->beg);
	free(b);
	sam_itr_destroy(data->iter);
return read_count;
}

// This function reads the regions file - generates the read count, depth and baf for each position
void get_cnv_data_for_regions_file (fdata **data, char *fname) {
  	char line[100];
  	float depth=0;
  	int i, j=1, read_count=0;
  	FILE *infile=fopen(fname, "r");
  	char *tl, *baf_bait="", *rd_bait="", *ref_base="", seps[] = "\t";
  	while(fgets(line, sizeof(line), infile)!=NULL) {
  		for (i = 0; i < strlen(line); i++) {
            if ( line[i] == '\n' || line[i] == '\r' )
                line[i] = '\0';
        }
    	for (tl = strtok (line, seps); tl; tl = strtok (NULL, seps))  {
      		if (j==1) baf_bait = tl;
      		if (j==2) ref_base = tl;
      		if (j==3) rd_bait = tl;
      		j++;
    	}j=1;
		float baf = extract_allele_counts(data, baf_bait, ref_base);
		read_count = read_cnv_data(data[0], rd_bait, &depth);
		char *gt = convert_baf_to_genotype(baf);
  		printf("%s\t%d\t%f\t%s\t%s\t%f\t%s\n", rd_bait, read_count, depth, baf_bait, ref_base, baf, gt);
  	}
 	fclose(infile);
}

int main_cnv(int argc, char *argv[]) {
	int n, mapQ=10;
	char *fname="regions_file";
    while ((n = getopt(argc, argv, "Q:R:E:")) >= 0) {
        switch (n) {
            case 'Q': mapQ = atoi(optarg); break;    // mapping quality threshold
            case 'R': fname = optarg; break;    // regions file
        }
    }
    if (optind == argc) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage: cnvdata [options] in1.bam or in1.cram\n");
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "   -Q <int>            mapping quality threshold\n");
        fprintf(stderr, "   -R <char*>          regions file\n");
        fprintf(stderr, "\n");
        return 1;
    }
	fdata **data;
    data = calloc(1, sizeof(fdata*));
    data[0] =get_fdata(argv[optind], mapQ);
	get_cnv_data_for_regions_file(data, fname);
    hts_idx_destroy(data[0]->idx);
    bam_hdr_destroy(data[0]->hdr);
    if (data[0]->fp) sam_close(data[0]->fp);
    // hts_itr_destroy(data[0]->iter);
    free(data);
return 0;
}

int main(int argc, char *argv[]) {
return main_cnv(argc, argv);
}
