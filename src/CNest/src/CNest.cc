#include <map>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cstring>
#include <math.h>
#include <fstream>
#include <iostream>

using namespace std;

#include <R.h> 
#include <Rdefines.h> 
#include <Rinternals.h>
#include <Rinterface.h>

/*
	Description:	CNest - the copy number estimator
	Authors:		Tomas Fitzgerald
	Date:			5/09/2011
*/

extern "C" {

// Definition of the depth_seq class - default construtors, values contained in class and any friends  
class depth_seq {
	public:
    	depth_seq(): d_count(0), d_depth(0), d_depth_ref(0), d_depth_err(0) {}
    	depth_seq(double count, double depth, double depth_ref, double depth_err): d_count(count), d_depth(depth), d_depth_ref(depth_ref), d_depth_err(depth_err){}
 
    friend ostream& operator<<(ostream&, const depth_seq&);
	
	public:
   		double d_count;
   		double d_depth;
   		double d_depth_ref;
   		double d_depth_err;
};

class ratio_seq {
    public:
        ratio_seq(): d_nratio(0), d_nratio_gen_matched(0), d_nratio_gen_mismatched(0) {}
        ratio_seq(double nratio, double nratio_gen_matched, double nratio_gen_mismatched): d_nratio(nratio), d_nratio_gen_matched(nratio_gen_matched), d_nratio_gen_mismatched(nratio_gen_mismatched){}

    public:
        double d_nratio;
        double d_nratio_gen_matched;
        double d_nratio_gen_mismatched;
};
    
// friend of depth_seq class - a ostream call on an depth_seq class to stout
ostream& operator<<(ostream& out, const depth_seq& c) {
   	printf("%f\t%f\t%f\t%f", c.d_count, c.d_depth, c.d_depth_ref, c.d_depth_err);
    return out;
}

int get_size_for_type(int item_type) {
	int item_size = 0;
	if(item_type==1) { depth_seq a(0,0,0,0); item_size=sizeof(a); }
	if(item_type==2) { ratio_seq a(0,0,0); item_size=sizeof(a); }
	return item_size;
}

// Read the index file and put the position values in our map - keyed by chr
map < int, vector<string> > readIndexFile(const char *fname) {
	map < int, vector<string> > index_data;
  	char line[100]; FILE *infile; infile=fopen(fname, "r"); 
  		while(fgets(line, sizeof(line), infile)!=NULL) {
  			line[strlen(line)-1]='\0';
    		string ss(line); string s=strtok(line, "\t");
    		istringstream buffer(s); int c; buffer >> c;
    		index_data[c].push_back(ss);
  		}
 	fclose(infile);
 	return index_data;
}

// Loop the index file - calcuate the starting byte position and the number of rows (classes) to read from the binary file
// Note: index file must be numerically sorted by chr, start and stop - (a map in c++ automatically sorts its keys numerically) 
vector < vector<int> > poffset(const char *fname, int chr, int start, int stop, int item_size) {
	map < int, vector<string> > index_data = readIndexFile(fname);
	long offset=0; int nrow=0; vector < vector<int> > positions;
	for(map<int, vector<string> >::const_iterator it = index_data.begin(); it != index_data.end(); it++) {
		int chrom=it->first; vector<string> features=it->second; int len=features.size();
		if(chrom != chr) {
			offset+=item_size*len;
		} else {
			for(int i=0;i<len;i++) {
				int count=0; int sta; int sto;
				string fea(features[i]); char* f=new char[fea.size()+1]; f[fea.size()]=0;
				memcpy(f,fea.c_str(),fea.size()); string ss=strtok(f, "\t");	
				while (f!=NULL) {
			 		int t=atoi(f);
			 		if(count ==1) { sta = t; }
			 		if(count ==2) { sto = t; }
    				f=strtok(NULL, "\t");
    				count++;
  				}					
  				if(sta < start) {
  					offset+=item_size;
  				} else if (sto <= stop) {
  					vector<int> v; v.push_back(chr); v.push_back(sta); v.push_back(sto);
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

// Utility function to create binary format from depth_seq input file format 
void depth_seq_bin(char** inputname, char** binname) {	
	string s; ifstream file( *inputname ); ofstream binfile(*binname, std::ios::out | std::ios::binary);
	while (std::getline(file, s)){
		vector<string> v; stringstream ss; ss << s;
  		while (getline( ss, s, '\t' )) v.push_back( s );
 		istringstream buffer1(v[3]); double ratio; buffer1 >> ratio;
 		istringstream buffer2(v[4]); double scan_s; buffer2 >> scan_s;
 		istringstream buffer3(v[5]); double red; buffer3 >> red;
 		istringstream buffer4(v[6]); double green; buffer4 >> green;
		depth_seq c(ratio, red, green, scan_s); binfile.write((char *)(&c), sizeof(c));
	}
	file.close(); binfile.close();
}
    
// Utility function to create binary format from ratio_seq input file format
void ratio_seq_bin(char** inputname, char** binname) {
        string s; ifstream file( *inputname ); ofstream binfile(*binname, std::ios::out | std::ios::binary);
        while (std::getline(file, s)){
            vector<string> v; stringstream ss; ss << s;
            while (getline( ss, s, '\t' )) v.push_back( s );
            istringstream buffer1(v[3]); double ratio; buffer1 >> ratio;
            istringstream buffer2(v[4]); double ratio_matched; buffer2 >> ratio_matched;
            istringstream buffer3(v[5]); double ratio_mismatched; buffer3 >> ratio_mismatched;
            ratio_seq c(ratio, ratio_matched, ratio_mismatched); binfile.write((char *)(&c), sizeof(c));
        }
        file.close(); binfile.close();
}

void get_offset(char** filename, int *chr, int *start, int *stop, int *rows, int *offset, int *item_type) {
	int item_size = get_size_for_type((int)*item_type);
	cout << item_size << endl;
	vector < vector<int> > positions = poffset(*filename, (int)*chr, (int)*start, (int)*stop, item_size);
	vector<int> v = positions[positions.size()-1];
	*offset = v[0]; *rows=v[1];
}

void get_positions(char** filename, int *chr, int *start, int *stop, int *chr_values, int *start_values, int *stop_values, int *item_type) {
	vector < vector<int> > positions = poffset(*filename, (int)*chr, (int)*start, (int)*stop, get_size_for_type((int)*item_type));
	vector<int> v = positions[positions.size()-1];
	for(int i=0;i<v[1];i++) {
 		vector<int> v = positions[i];
 		chr_values[i] = v[0]; start_values[i] = v[1]; stop_values[i] = v[2];
 	}
}

void get_values(char** filename, int *starting_pos, int *number_row_to_read, double *values, int *item_type, int *return_type) {
	long s_pos = (long)*starting_pos;
	cout << s_pos << endl;
	cout << (int)*item_type << endl;
	ifstream file(*filename, ios::binary); file.seekg(s_pos);
	if((int)*item_type == 1) {
		depth_seq c(0,0,0,0);
		for(int i=0;i<(int)*number_row_to_read;i++) {
			file.read((char *)(&c), sizeof(c)); 
			if((int)*return_type == 1) { values[i] = c.d_count; }
			if((int)*return_type == 2) { values[i] = c.d_depth; }
			if((int)*return_type == 3) { values[i] = c.d_depth_ref; }
			if((int)*return_type == 4) { values[i] = c.d_depth_err; }
		}
	}
	if((int)*item_type == 2) {
		ratio_seq c(0,0,0);
		for(int i=0;i<(int)*number_row_to_read;i++) {
			file.read((char *)(&c), sizeof(c)); 
			if((int)*return_type == 1) { values[i] = c.d_nratio; }
			if((int)*return_type == 2) { values[i] = c.d_nratio_gen_matched; }
			if((int)*return_type == 3) { values[i] = c.d_nratio_gen_mismatched; }
		}
	}
	file.close(); 
}	

void compress_rbin(char** binname, char** filenames, int *n, int *starting_pos, int *number_row_to_read) {
    ratio_seq c; 
    ofstream binfile(*binname, std::ios::out | std::ios::binary);
    long start_pos = (long)*starting_pos;
    for(int i=0;i<(int)*number_row_to_read;i++) {
    	long toff = sizeof(c) * i;
	    for(int j=0;j<(int)*n;j++) {
		    char* filename = filenames[j];
		    ifstream file(filename, ios::binary); 
		    file.seekg(start_pos+toff);
		 	file.read((char *)(&c), sizeof(c)); 
		 	binfile.write((char *)(&c), sizeof(c));
		 	file.close(); 
		 }
	 }
	 binfile.close();
}

void write_blank_cbin_file(char** binname, int *n, int *number_row_to_read) {
    ratio_seq c(0, 0, 0);
    ofstream binfile(*binname, std::ios::out | std::ios::binary);
    for(int i=0;i<(int)*number_row_to_read;i++) {
	    for(int j=0;j<(int)*n;j++) {
		 	binfile.write((char *)(&c), sizeof(c));
		 }
	 }
	 binfile.close();
}

void exome_ratio_read_write_chunk_array(char** binname, char** filenames, int *n_values, int *chunk_size, int *n_samples, int *starting_pos) {
	int sample_len = (int)*n_samples;
	int n_loops = ceil( (int)*n_values / (int)*chunk_size )+1;
	cout << n_loops << endl;
	long start_pos = (long)*starting_pos;
	ofstream binfile(*binname, std::ios::out | std::ios::binary);
	for(int k=1;k<=n_loops;k++) {
		int new_vals = *chunk_size;
		cout << k << endl;
		if(((int)*chunk_size)*k > (int)*n_values) { new_vals= (int)*n_values - ((int)*chunk_size)*(k-1); }
		ratio_seq c; 
		int nlen = (new_vals*sample_len);
		ratio_seq cs_out[nlen];
		for(int j=0;j<sample_len;j++) {
			char* filename = filenames[j];
			ifstream file(filename, ios::binary); 
			file.seekg(start_pos);
			for(int i=0;i<new_vals;i++) {
				file.read((char *)(&c), sizeof(c)); 
				int pin = j+(i*sample_len); 
				cs_out[pin] = c;
			}
			file.close(); 
		}
    	binfile.write((char *)(&cs_out), sizeof(cs_out));
		long toff = sizeof(c) * ((int)*chunk_size); 
		start_pos = start_pos+toff;
	}
	binfile.close();
}

void write_sample_to_blank_cbin(char** binname, char** filename, int *sample_index, int *number_samples, int *starting_pos, int *number_row_to_read) {
	ratio_seq c;
    ifstream file(*filename, ios::binary); 
    fstream binfile(*binname, std::ios::binary|std::ios::out|std::ios::in);
   	long bin_pos = ((*sample_index)*sizeof(c))-sizeof(c);
	file.seekg((long)*starting_pos);
	for(int i=0;i<(int)*number_row_to_read;i++) {
		long toff = ((*number_samples) * sizeof(c)) * i;
		binfile.seekg(bin_pos+toff);
		file.read((char *)(&c), sizeof(c)); 
		binfile.write((char *)(&c), sizeof(c));
	}
	file.close(); 
	binfile.close();
}

void get_cbin_index_pos(char** filename, int *index_pos, int *number_samples, int *cntype, double *values) {
    ratio_seq c;
    long starting_pos = sizeof(c) * (long)((*index_pos-1)*(*number_samples));
    ifstream file(*filename, ios::binary); file.seekg(starting_pos);
    for(int i=0;i<(int)*number_samples;i++) {
 		file.read((char *)(&c), sizeof(c)); 
 		if(*cntype==1) { values[i] = c.d_nratio; }
 		if(*cntype==2) { values[i] = c.d_nratio_gen_matched; }
 		if(*cntype==3) { values[i] = c.d_nratio_gen_mismatched; }
 	}
    file.close(); 
}


void get_cbin_index_range(char** filename, int *index_pos, int *number_samples, int *number_positions, int *cntype, double *values) {
    ratio_seq c;
    long starting_pos = sizeof(c) * (long)((*index_pos-1)*(*number_samples));
    ifstream file(*filename, ios::binary); file.seekg(starting_pos); int pin=0;
    for(int j=0;j<(int)*number_positions;j++) {
	    for(int i=0;i<(int)*number_samples;i++) {
	 		file.read((char *)(&c), sizeof(c)); 
	 		if(*cntype==1) { values[pin] = c.d_nratio; }
	 		if(*cntype==2) { values[pin] = c.d_nratio_gen_matched; }
	 		if(*cntype==3) { values[pin] = c.d_nratio_gen_mismatched; }
	 		pin++;
	 	}
	}
    file.close(); 
}

	
}
