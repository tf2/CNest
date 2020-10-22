#ifndef CNVDATA_H
#define CNVDATA_H

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include "/htslib/htslib/sam.h" 
#include "samtools.h"

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

fdata *get_fdata (char *filename, int mapQ);
int read_cnv_data (fdata *data, char *reg, float *depth) ;;
void write_cnv_data_for_regions_file (fdata *data, char *fname, char *bname);
int count_lines(char *fname);
bdata** get_cnv_data_for_regions_file (fdata *data, char *fname);
//bdata** get_bdata (char *bname, int starting_pos, int number_row_to_read);

#ifdef __cplusplus
} // closing brace for extern "C"

#endif
#endif
