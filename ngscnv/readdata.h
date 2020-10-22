#ifndef READDATA_H
#define READDATA_H

#include <map>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cstring>
#include <math.h>
#include <fstream>
#include <iostream>

using namespace std;


class readdata {
	public:

    	readdata(): d_readcount(0), d_depth(0), d_ref_readcount(0), d_ref_mad(0), d_l2r(0), d_nl2r(0) {}
    	readdata(double readcount, double depth): d_readcount(readcount), d_depth(depth), d_ref_readcount(0), d_ref_mad(0), d_l2r(0), d_nl2r(0) {}
    	readdata(double readcount, float depth, int ref_readcount, double ref_mad, double l2r, double nl2r): d_readcount(readcount), d_depth(depth), d_ref_readcount(ref_readcount), d_ref_mad(ref_mad), d_l2r(l2r), d_nl2r(nl2r){}
 
	public:
		int d_readcount;
    	double d_depth;
    	int d_ref_readcount;
    	double d_ref_mad;
    	double d_l2r;
    	double d_nl2r;
};

map < int, vector<string> > readIndexFile(const char *fname);
vector < vector<int> > findeOffset(const char *fname, int chr, int start, int stop);
void write_bam_to_rd(char *bamname, char *fname, char *bname);
vector < vector<int> > findOffset(const char *fname, string chr, int start, int stop);
vector < int > getOffsets(char* filename, string chr, int start, int stop);
vector < readdata > getValues(char* filename, int starting_pos, int number_row_to_read);

#endif