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
	Description:	Binary Writing and Reading Methods for Fast Data Access from R
	Authors:		Tomas Fitzgerald
	Date:			5/09/2011
*/

extern "C" {

//int OFFSET; int NUMROW;
//map < int, vector<string> > index_data;

// Definition of the aCGH class - default construtors, values contained in class and any friends  
class aCGH {
	public:
    	aCGH(): d_ratio(0), d_red(0), d_green(0), d_scan_s(0) {}
    	aCGH(double ratio, double red, double green, double scan_s): d_ratio(ratio), d_red(red), d_green(green), d_scan_s(scan_s){}
 
    friend ostream& operator<<(ostream&, const aCGH&);
	
	public:
   		double d_ratio;
   		double d_red;
   		double d_green;
   		double d_scan_s;
};

class exome {
	public:
    	exome(): d_nratio(0), d_adm3score(0) {}
    	exome(double nratio, double adm3score): d_nratio(nratio), d_adm3score(adm3score){}
 
    //friend ostream& operator<<(ostream&, const exome&);
	
	public:
		double d_nratio;
		double d_adm3score;
};

class exome_ratio_reference {
    public:
        exome_ratio_reference(): d_nratio(0), d_nratio_gen_matched(0), d_nratio_gen_mismatched(0) {}
        exome_ratio_reference(double nratio, double nratio_gen_matched, double nratio_gen_mismatched): d_nratio(nratio), d_nratio_gen_matched(nratio_gen_matched), d_nratio_gen_mismatched(nratio_gen_mismatched){}
        
        //friend ostream& operator<<(ostream&, const exome&);
        
    public:
        double d_nratio;
        double d_nratio_gen_matched;
        double d_nratio_gen_mismatched;
};
    
    
// friend of aCGH class - when we make a ostream call on an aCGH class (e.g. using cout) - we will print it in the following way 
ostream& operator<<(ostream& out, const aCGH& c) {
   	printf("%f\t%f\t%f\t%f", c.d_ratio, c.d_red, c.d_green, c.d_scan_s);
    return out;
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
vector < vector<int> > findOffset(const char *fname, int chr, int start, int stop) {
	map < int, vector<string> > index_data = readIndexFile(fname);
	aCGH a; long offset=0; int nrow=0; vector < vector<int> > positions;
	for(map<int, vector<string> >::const_iterator it = index_data.begin(); it != index_data.end(); it++) {
		int chrom=it->first; vector<string> features=it->second; int len=features.size();
		if(chrom != chr) {
			offset+=sizeof(a)*len;
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
  					offset+=sizeof(a);
  				} else if (sto <= stop) {
  					vector<int> v; v.push_back(chr); v.push_back(sta); v.push_back(sto);
  					positions.push_back(v); nrow++;
					//offset+=sizeof(a);
  					if(sto>stop) { break; }
  				} else {
  					break;
  				}
			}
			vector<int> v; v.push_back(offset); v.push_back(nrow);
  			positions.push_back(v);
			//OFFSET = offset; NUMROW=nrow;
			break;
		}
	}

	return positions;
}

// Utility function to create binary format from aCGH input file format 
void makeaCGHbin(char** inputname, char** binname) {	
	string s; double ratio, red, green, scan_s;
	ifstream file( *inputname ); ofstream binfile(*binname, std::ios::out | std::ios::binary);
	while (std::getline(file, s)){
		vector<string> v; stringstream ss; ss << s;
  		while (getline( ss, s, '\t' )) v.push_back( s );
 		istringstream buffer1(v[3]); double ratio; buffer1 >> ratio;
 		istringstream buffer2(v[4]); double scan_s; buffer2 >> scan_s;
 		istringstream buffer3(v[5]); double red; buffer3 >> red;
 		istringstream buffer4(v[6]); double green; buffer4 >> green;
		aCGH c(ratio, red, green, scan_s); binfile.write((char *)(&c), sizeof(c));
	}
	file.close(); binfile.close();
}
    
// Utility function to create binary format from aCGH input file format
void make_exome_ratio_reference_bin(char** inputname, char** binname) {
        string s; double ratio, ratio_matched, ratio_mismatched;
        ifstream file( *inputname ); ofstream binfile(*binname, std::ios::out | std::ios::binary);
        while (std::getline(file, s)){
            vector<string> v; stringstream ss; ss << s;
            while (getline( ss, s, '\t' )) v.push_back( s );
            istringstream buffer1(v[3]); double ratio; buffer1 >> ratio;
            istringstream buffer2(v[4]); double scan_s; buffer2 >> ratio_matched;
            istringstream buffer3(v[5]); double red; buffer3 >> ratio_mismatched;
            exome_ratio_reference c(ratio, ratio_matched, ratio_mismatched); binfile.write((char *)(&c), sizeof(c));
        }
        file.close(); binfile.close();
}
    
    
	// Loop the index file - calcuate the starting byte position and the number of rows (classes) to read from the binary file
	// Note: index file must be numerically sorted by chr, start and stop - (a map in c++ automatically sorts its keys numerically) 
	vector < vector<int> > findeOffset(const char *fname, int chr, int start, int stop) {
		map < int, vector<string> > index_data = readIndexFile(fname);
		exome a; long offset=0; int nrow=0; vector < vector<int> > positions;
		for(map<int, vector<string> >::const_iterator it = index_data.begin(); it != index_data.end(); it++) {
			int chrom=it->first; vector<string> features=it->second; int len=features.size();
			if(chrom != chr) {
				offset+=sizeof(a)*len;
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
  					delete f;
					if(sta < start) {
						offset+=sizeof(a);
					} else if (sto <= stop) {
						vector<int> v; v.push_back(chr); v.push_back(sta); v.push_back(sto);
						positions.push_back(v); nrow++;
						//offset+=sizeof(a);
						if(sto>stop) { break; }
					} else {
						break;
					}
				}
				vector<int> v; v.push_back(offset); v.push_back(nrow);
				positions.push_back(v);
				//OFFSET = offset; NUMROW=nrow;
				break;
			}
		}
		
		return positions;
	}

    
	// Loop the index file - calcuate the starting byte position and the number of rows (classes) to read from the binary file
	// Note: index file must be numerically sorted by chr, start and stop - (a map in c++ automatically sorts its keys numerically) 
vector < vector<int> > findeOffset_ratio(const char *fname, int chr, int start, int stop) {
		map < int, vector<string> > index_data = readIndexFile(fname);
		exome_ratio_reference a; long offset=0; int nrow=0; vector < vector<int> > positions;
		for(map<int, vector<string> >::const_iterator it = index_data.begin(); it != index_data.end(); it++) {
			int chrom=it->first; vector<string> features=it->second; int len=features.size();
			if(chrom != chr) {
				offset+=sizeof(a)*len;
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
  					delete f;
					if(sta < start) {
						offset+=sizeof(a);
					} else if (sto <= stop) {
						vector<int> v; v.push_back(chr); v.push_back(sta); v.push_back(sto);
						positions.push_back(v); nrow++;
						//offset+=sizeof(a);
						if(sto>stop) { break; }
					} else {
						break;
					}
				}
				vector<int> v; v.push_back(offset); v.push_back(nrow);
				positions.push_back(v);
				//OFFSET = offset; NUMROW=nrow;
				break;
			}
		}
		
		return positions;
}

	
	
// Utility function to create binary format from exome input file format 
void makeexomebin(char** inputname, char** binname) {	
	string s; double nratio, ratio, error, depth;
	ifstream file( *inputname ); ofstream binfile(*binname, std::ios::out | std::ios::binary);
	while (std::getline(file, s)){
		vector<string> v; stringstream ss; ss << s;
  		while (getline( ss, s, '\t' )) v.push_back( s );
 		istringstream buffer1(v[5]); double nratio; buffer1 >> nratio;
 		istringstream buffer2(v[7]); double adm3score; buffer2 >> adm3score;
		exome c(nratio, adm3score); binfile.write((char *)(&c), sizeof(c));
	}
	file.close(); binfile.close();
}

// First function to call from R 
// - it will search the index file returning the length of data to extract and the starting byte position
// - we only need do this once for muliple file reads
void getOffsets(char** filename, int *chr, int *start, int *stop, int *rows, int *offset) {
	//readIndexFile(*filename); 
	vector < vector<int> > positions = findOffset(*filename, (int)*chr, (int)*start, (int)*stop);
	// Here we set the pointers from R to the values we need
	vector<int> v = positions[positions.size()-1];
	*offset = v[0]; *rows=v[1];
	//*rows=NUMROW; *offset = OFFSET;
}
	
// First function to call from R 
// - it will search the index file returning the length of data to extract and the starting byte position
// - we only need do this once for muliple file reads
void geteOffsets(char** filename, int *chr, int *start, int *stop, int *rows, int *offset) {
	//readIndexFile(*filename); 
	vector < vector<int> > positions = findeOffset(*filename, (int)*chr, (int)*start, (int)*stop);
	// Here we set the pointers from R to the values we need
	vector<int> v = positions[positions.size()-1];
	*offset = v[0]; *rows=v[1];
	//*rows=NUMROW; *offset = OFFSET;
}


// First function to call from R 
// - it will search the index file returning the length of data to extract and the starting byte position
// - we only need do this once for muliple file reads
void geteOffsets_ratio(char** filename, int *chr, int *start, int *stop, int *rows, int *offset) {
	//readIndexFile(*filename); 
	vector < vector<int> > positions = findeOffset_ratio(*filename, (int)*chr, (int)*start, (int)*stop);
	// Here we set the pointers from R to the values we need
	vector<int> v = positions[positions.size()-1];
	*offset = v[0]; *rows=v[1];
	//*rows=NUMROW; *offset = OFFSET;
}

void getePositions_ratio(char** filename, int *chr, int *start, int *stop, int *chr_values, int *start_values, int *stop_values) {
	//readIndexFile(*filename); 
	vector < vector<int> > positions = findeOffset_ratio(*filename, (int)*chr, (int)*start, (int)*stop);
	vector<int> v = positions[positions.size()-1];
	for(int i=0;i<v[1];i++) {
 		vector<int> v = positions[i];
 		chr_values[i] = v[0]; start_values[i] = v[1]; stop_values[i] = v[2];
 	}
}

void geteValues1_ratio(char** filename, int *starting_pos, int *number_row_to_read, double *values) {
    exome_ratio_reference c; ifstream file(*filename, ios::binary); file.seekg((long)*starting_pos);
    for(int i=0;i<(int)*number_row_to_read;i++) {
 		file.read((char *)(&c), sizeof(c)); 
 		values[i] = c.d_nratio;
 	}
    file.close(); 
}

void geteValues2_ratio(char** filename, int *starting_pos, int *number_row_to_read, double *values) {
	exome_ratio_reference c; ifstream file(*filename, ios::binary); file.seekg((long)*starting_pos);
	for(int i=0;i<(int)*number_row_to_read;i++) {
		file.read((char *)(&c), sizeof(c)); 
		values[i] = c.d_nratio_gen_matched;
	}
	file.close(); 
}
	

void geteValues3_ratio(char** filename, int *starting_pos, int *number_row_to_read, double *values) {
	exome_ratio_reference c; ifstream file(*filename, ios::binary); file.seekg((long)*starting_pos);
	for(int i=0;i<(int)*number_row_to_read;i++) {
		file.read((char *)(&c), sizeof(c)); 
		values[i] = c.d_nratio_gen_mismatched;
	}
	file.close(); 
}
	




// Second function to call from R 
// - here we additionally parse in pointers to new integer arrays of the correct size (NUMROW) from R 
// - we fill these arrays with the chr, start and stop positions from the index file
// - again we only need do this once.
void getPositions(char** filename, int *chr, int *start, int *stop, int *chr_values, int *start_values, int *stop_values) {
	//readIndexFile(*filename); 
	vector < vector<int> > positions = findOffset(*filename, (int)*chr, (int)*start, (int)*stop);
	vector<int> v = positions[positions.size()-1];
	for(int i=0;i<v[1];i++) {
 		vector<int> v = positions[i];
 		chr_values[i] = v[0]; start_values[i] = v[1]; stop_values[i] = v[2];
 	}
}

// Thrid function to call from R 
// - we read the binary file given the starting byte position and number of rows, filling a double array with the values 
// - we can do this as many times as we like (i.e. for different DDD samples returning a new double array to R each time)
void getValues(char** filename, int *starting_pos, int *number_row_to_read, double *values) {
	aCGH c; ifstream file(*filename, ios::binary); file.seekg((long)*starting_pos);
	for(int i=0;i<(int)*number_row_to_read;i++) {
		file.read((char *)(&c), sizeof(c)); 
		values[i] = c.d_ratio;
	}
	file.close(); 
}	
	
void geteValues1(char** filename, int *starting_pos, int *number_row_to_read, double *values) {
    exome c; ifstream file(*filename, ios::binary); file.seekg((long)*starting_pos);
    for(int i=0;i<(int)*number_row_to_read;i++) {
 		file.read((char *)(&c), sizeof(c)); 
 		values[i] = c.d_nratio;
 	}
    file.close(); 
}

void geteValues2(char** filename, int *starting_pos, int *number_row_to_read, double *values) {
	exome c; ifstream file(*filename, ios::binary); file.seekg((long)*starting_pos);
	for(int i=0;i<(int)*number_row_to_read;i++) {
		file.read((char *)(&c), sizeof(c)); 
		values[i] = c.d_adm3score;
	}
	file.close(); 
}

void compress_rbin(char** binname, char** filenames, int *n, int *starting_pos, int *number_row_to_read) {
    exome_ratio_reference c; 
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
    exome_ratio_reference c(0, 0, 0);
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
		exome_ratio_reference c; 
		int nlen = (new_vals*sample_len);
		exome_ratio_reference cs_out[nlen];
		for(int j=0;j<sample_len;j++) {
			char* filename = filenames[j];
			cout << j << endl;
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
	exome_ratio_reference c;
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
    exome_ratio_reference c;
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

	
}
