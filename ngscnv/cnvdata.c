/* This program to extract copy number data from BAM / CRAM files
 * Author: Tomas William Fitzgerald
 * To compile this program, you may:
 *
 *   gcc -O3 -Wall -Isrc -I/nfs/ddd0/Tom/cram/cram/samtools-1.2/ -I/nfs/ddd0/Tom/cram/cram/samtools-1.2/htslib-1.2.1/ -I/nfs/ddd0/Tom/cram/cram/samtools-1.2/htslib-1.2.1/htslib -rdynamic -o cnvdata -L/nfs/ddd0/Tom/cram/cram/samtools-1.2/htslib-1.2.1/ -L/nfs/ddd0/Tom/cram/cram/samtools-1.2/ cnvdata.c -lhts -lbam -lpthread -lz -lm

 */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include "/htslib/htslib/sam.h" 
#include "samtools.h"
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>

typedef struct {     // file parsing information
    htsFile *fp;     // the file handle
    bam_hdr_t *hdr;  // the file header
    hts_itr_t *iter; // NULL if a region not specified
    hts_idx_t *idx; // the hts index
    int min_mapQ;
} fdata;

typedef struct {     // data storing info
    int readcount;
    float depth;
    int ref_readcount;
    double ref_mad;
    double l2r;
    double nl2r;
} bdata;

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

int read_cnv_data (fdata *data, char *reg, float *depth) {
	int i, result, read_count=0, coverage = 0;
	bam1_t* b = NULL; b = bam_init1();
	data->iter = sam_itr_querys(data->idx, data->hdr, reg); // set the iterator	
	while ((result = sam_itr_next(data->fp, data->iter, b)) >= 0) {
    	/*if(b->core.qual < data->min_mapQ || (b->core.flag & BAM_FUNMAP)
			|| !(b->core.flag & BAM_FPROPER_PAIR) || (b->core.flag & BAM_FMUNMAP)//Proper pair and mate unmapped
			|| (b->core.flag & BAM_FDUP)//1024 is PCR/optical duplicate
			|| (b->core.flag & BAM_FSECONDARY) || (b->core.flag & BAM_FQCFAIL)//Secondary alignment, quality fail
			|| (b->core.flag & BAM_FSUPPLEMENTARY) ) continue; */
    	for(i=0;i<b->core.l_qseq;i++) {
    		if(b->core.pos >= data->iter->beg && b->core.pos+i <= data->iter->end)
    			coverage += 1;
    	}
    	read_count++;
	}
	*depth = (float)coverage /(data->iter->end-data->iter->beg);
return read_count;
}

void write_cnv_data_for_regions_file (fdata *data, char *fname, char *bname) {
	bdata readdata;
  	char line[100];
  	float depth=0;
  	int i, read_count=0;
  	FILE *infile=fopen(fname, "rb");
  	FILE *fp;
	fp=fopen(bname, "w");
  	while(fgets(line, sizeof(line), infile)!=NULL) {
  		for (i = 0; i < strlen(line); i++) {
            if ( line[i] == '\n' || line[i] == '\r' )
                line[i] = '\0';
        }
		read_count = read_cnv_data(data, strdup(line), &depth);
  		printf("%s\t%d\t%f\n", line, read_count, depth);
  		readdata.readcount = read_count;
  		readdata.depth = depth;
  		
  		fwrite(&readdata, sizeof(readdata), 1, fp);
  	}
 	fclose(infile);
 	fclose(fp);
}

int count_lines(char *fname) {
	int count=0;
	char line[100];
	FILE *infile=fopen(fname, "r");
  	while(fgets(line, sizeof(line), infile)!=NULL)
		count++;
 	fclose(infile);
return count;
}

bdata** get_cnv_data_for_regions_file (fdata *data, char *fname) {
	int number_lines = count_lines(fname);
	bdata** readdatas = malloc(number_lines * sizeof(*readdatas));
  	char line[100];
  	float depth=0;
  	int i, k=0, read_count=0;
  	FILE *infile=fopen(fname, "rb");
  	while(fgets(line, sizeof(line), infile)!=NULL) {
  		for (i = 0; i < strlen(line); i++) {
            if ( line[i] == '\n' || line[i] == '\r' )
                line[i] = '\0';
        }
		read_count = read_cnv_data(data, strdup(line), &depth);
		readdatas[k] = malloc(sizeof(**readdatas));
		readdatas[k]->readcount = read_count;
		readdatas[k]->depth = depth;
		k++;
  	}
 	fclose(infile);
return readdatas;
}

