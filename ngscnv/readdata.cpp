#include <map>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cstring>
#include <math.h>
#include <fstream>
#include <iostream>

using namespace std;

extern "C" {
	#include "cnvdata.h"
}

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

void write_bam_to_rd(char *bamname, char *fname, char *bname) {
	string line;
	float depth=0;
	ifstream file(fname); 
	ofstream binfile(bname, std::ios::out | std::ios::binary);
	fdata *data = get_fdata(bamname, 0);
	while (std::getline(file, line)){
		int read_count = read_cnv_data(data, strdup(line.c_str()), &depth);
		readdata c(read_count, depth);
		binfile.write((char *)(&c), sizeof(c));
		//cout << line << "\t" << read_count << "\t" << depth << endl;
	}
	file.close(); binfile.close();
}

// Read the index file and put the position values in our map - keyed by chr
map < string, vector<string> > readIndexFile(const char *fname) {
	map < string, vector<string> > index_data;
  	char line[100]; FILE *infile; infile=fopen(fname, "r"); 
  		while(fgets(line, sizeof(line), infile)!=NULL) {
  			line[strlen(line)-1]='\0';
  			stringstream test(line);
			string segment;
			vector<string> seglist;
			while(std::getline(test, segment, ':')) {
   				seglist.push_back(segment);
			}

    		index_data[seglist[0]].push_back(seglist[1]);
  		}
 	fclose(infile);
 	return index_data;
}

// Loop the index file - calcuate the starting byte position and the number of rows (classes) to read from the binary file
// Note: index file must be numerically sorted by chr, start and stop - (a map in c++ automatically sorts its keys numerically) 
vector < vector<int> > findOffset(const char *fname, string chr, int start, int stop) {
	map < string, vector<string> > index_data = readIndexFile(fname);
	readdata a; long offset=0; int nrow=0; vector < vector<int> > positions;
	for(map<string, vector<string> >::const_iterator it = index_data.begin(); it != index_data.end(); it++) {
		string chrom=it->first; vector<string> features=it->second; int len=features.size();
		if(chrom != chr) {
			offset+=sizeof(a)*len;
		} else {
			for(int i=0;i<len;i++) {
				int count=0; int sta; int sto;
				string fea(features[i]); char* f=new char[fea.size()+1]; f[fea.size()]=0;
				memcpy(f,fea.c_str(),fea.size()); string ss=strtok(f, "-");	
				while (f!=NULL) {
					int t=atoi(f);
					if(count ==0) { sta = t; }
					if(count ==1) { sto = t; }
    				f=strtok(NULL, "\t");
    				count++;
  				}	
  				delete f;
				if(sta < start) {
					offset+=sizeof(a);
				} else if (sto <= stop) {
					vector<int> v; v.push_back(atoi(chr.c_str())); 
					v.push_back(sta); v.push_back(sto);
					positions.push_back(v); nrow++;
					if(sto>stop) { break; }
				} else {
					break;
				}
			}
			vector<int> v; v.push_back(offset); v.push_back(nrow);
			positions.push_back(v);
			break;
		}
	}
		
	return positions;
}

vector < int > getOffsets(char* filename, string chr, int start, int stop) {
	vector<int> offset_and_number_rows;
	vector < vector<int> > positions = findOffset(filename, chr, start, stop);
	vector<int> v = positions[positions.size()-1];
	offset_and_number_rows.push_back(v[0]); // offset 
	offset_and_number_rows.push_back(v[1]); // number rows 

return offset_and_number_rows;
}

vector < readdata > getValues(char* filename, int starting_pos, int number_row_to_read) {
    readdata c; 
    vector < readdata > readdata_values;
    cout << filename << endl;
    ifstream file(filename, ios::binary); 
    file.seekg((long)starting_pos);
    for(int i=0;i<number_row_to_read;i++) {
 		file.read((char *)(&c), sizeof(c));
 		readdata_values.push_back(c);
 	}
    file.close(); 
    
return readdata_values;
}



