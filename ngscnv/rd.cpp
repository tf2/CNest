#include <map>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cstring>
#include <math.h>
#include <fstream>
#include <iostream>
#include <cmath>  
#include <util.h>
#include <fh.h>
#include <gc.h>
#include <cwavef.h>

using namespace std;

class Comp {
	vector<double>& _v;
  	public:
    	Comp(vector<double>& v) : _v(v) {}
    	bool operator()(size_t i, size_t j){
        	return _v[i] > _v[j];
   		}
};

class rd {
	public:
    	rd(): d_readdepth(0) {}
    	rd(int readcount, double readdepth): d_readcount(readcount), d_readdepth(readdepth){}
 
    friend ostream& operator<<(ostream&, const rd&);
	
	public:
		int d_readcount;
		double d_readdepth;
		double d_log2ratio;
		double g_log2ratio;
		double w_log2ratio;
		double s_log2ratio;
		double d_admscore;
};

ostream& operator<<(ostream& out, const rd& c) {
   	printf("%f\t%f\t%f\t%f\t%f\n", c.d_readdepth, c.d_log2ratio, c.g_log2ratio, c.w_log2ratio, c.s_log2ratio, c.d_admscore);
    return out;
}

void makerdbin(vector<string> filenames, string binname) {	
	string s;
	ofstream binfile(binname.c_str(), std::ios::out | std::ios::binary);
	for(int i=0;i<filenames.size();i++) {
		std::ifstream myfile(filenames[i].c_str());
		while (std::getline(myfile, s)) {
			stringstream ss; ss << s;
  			vector<string> v;
  			while (getline( ss, s, '\t' )) v.push_back( s );
			istringstream buffer(v[3]);
 			int cvalue; buffer >> cvalue;
 			istringstream buffer1(v[4]);
 			double rvalue; buffer1 >> rvalue;
 			rd rd_member(cvalue, rvalue);
 			binfile.write((char *)(&rd_member), sizeof(rd_member));
		}
		//cout << i << endl;
	}
	binfile.close();
}

vector < rd > get_sample (string binname, vector<string> filenames, int filename_position, int number_lines) {
	rd c;
	vector< rd > sample;
	ifstream file(binname.c_str(), ios::binary); 
	long starting_position = filename_position*(sizeof(c)*number_lines);
	file.seekg(starting_position);
	for(int i=0; i<number_lines;i++) {
		long npos = starting_position+(sizeof(c)*i);
		file.seekg(npos);
		file.read((char *)(&c), sizeof(c)); 
		sample.push_back(c);
	}
	file.close();
return sample;
}

vector < rd > get_probe (string binname, vector<string> filenames, int probe_position, int number_lines) {
	rd c;
	vector< rd > sample;
	ifstream file(binname.c_str(), ios::binary); 
	long sample_length = sizeof(c)*number_lines;
	int starting_position = sizeof(c)*probe_position;
	file.seekg(starting_position);
	for(int i=0; i<filenames.size();i++) {
		long npos = starting_position+(sample_length*i);
		file.seekg(npos);
		file.read((char *)(&c), sizeof(c)); 
		sample.push_back(c);
	}
	file.close();
	return sample;
}

/*
vector < rd > get_mean_probes_ratio (string binname, vector<string> filenames, int probe_position, int number_lines, int range_start, int range_stop) {
	rd c;
	vector< double > sample_means;
	ifstream file(binname.c_str(), ios::binary); 
	long sample_length = sizeof(c)*number_lines;
	int starting_position = sizeof(c)*probe_position;
	file.seekg(starting_position);
	for(int i=0; i<filenames.size();i++) {
		long npos = starting_position+(sample_length*i);
		//file.seekg(npos);
		vector< double > dd;
		for(int k=range_start;k<range_stop;k++) {
			long tpos = npos+(sizeof(c)*k);
			file.seekg(tpos);
			file.read((char *)(&c), sizeof(c)); 
			dd.push_back(c.d_log2ratio);
		}
		sample_means.push_back( mean ( dd ) );
	}
	file.close();
	return sample_means;
}
*/

vector < double > get_sample_readdepth (string binname, vector<string> filenames, int filename_position, int number_lines) {
	rd c;
	vector< double > depths;
	ifstream file(binname.c_str(), ios::binary); 
	long starting_position = filename_position*(sizeof(c)*number_lines);
	file.seekg(starting_position);
	for(int i=0; i<number_lines;i++) {
		long npos = starting_position+(sizeof(c)*i);
		file.seekg(npos);
		file.read((char *)(&c), sizeof(c)); 
		depths.push_back(c.d_readdepth+0.01);
	}
	file.close();
return depths;
}

vector <string> get_highest_correlated_filenames (string binname, int number_lines, vector<string> filenames, int index, int max_c_filenames) {
	vector < string > fns; vector < double > cors;
	vector < double >  mdepths = get_sample_readdepth(binname, filenames, index, number_lines);
	for(int j=0;j<filenames.size();j++) {
		if (index==j) continue;
		vector < double >  depths = get_sample_readdepth(binname, filenames, j, number_lines);
		double cp = corpearson (mdepths, depths); 
		cors.push_back(cp);
		fns.push_back(filenames[j]);
		cout << j << endl;
	}
	vector<double> indices;
	for(int i=0;i<cors.size();i++) {
		indices.push_back(i);
	}
	vector < string > correlated_samples;
	std::sort(indices.begin(), indices.end(), Comp(cors));
	for(int i=0;i<indices.size();i++) {
		if(correlated_samples.size() < max_c_filenames) {
			correlated_samples.push_back(fns[indices[i]]);
			//cout << fns[indices[i]] << endl;
		}
	}
	
return correlated_samples;
}

vector<int> get_sample_indexes (vector<string> filenames, vector<string> searchnames) {
	vector<int> sample_indexes;
	for(int j=0;j<searchnames.size();j++) {
		for(int i=0;i<filenames.size();i++) {
			if(searchnames[j]==filenames[i]) {
				sample_indexes.push_back(i);
			}
		}
	}
return sample_indexes;
}

vector<double> get_median_readdepth (string binname, vector<string> filenames, vector<int> filename_positions, int number_lines) {
	rd c;
	ifstream file(binname.c_str(), ios::binary); 
	vector<double> median_dp;
	for(int i=0; i<number_lines;i++) {
		long position = sizeof(c)*i;
		vector<double> med_depths;
		for(int j=0; j<filename_positions.size();j++) {
			long this_position = position + filename_positions[j]*(sizeof(c)*number_lines);
			file.seekg(this_position);
			file.read((char *)(&c), sizeof(c)); 
			med_depths.push_back(c.d_readdepth+0.01);
			//cout << c.d_readdepth+0.01 << endl;
		}
		double med = median(med_depths);
		//double med = find_maxima(med_depths, floor(med_depths.size()/5));
		median_dp.push_back(med);
	}
	file.close();
return median_dp;
}

vector <double> update_log2_ratio (string binname, vector<string> filenames, vector<double> median_dp, int filename_position, int number_lines) {
	rd c;
	vector<double> lr2_values;
	ifstream file(binname.c_str(), ios::binary);
	ofstream ofile(binname.c_str(), ios::binary | ios::out | ios::in);
	long starting_position = filename_position*(sizeof(c)*number_lines);
	file.seekg(starting_position);
	for(int i=0; i<number_lines;i++) {
		long npos = starting_position+(sizeof(c)*i);
		file.seekg(npos);
		file.read((char *)(&c), sizeof(c));
		double l2 = log2 ( (c.d_readdepth+0.01)/median_dp[i] );
		//cout << c.d_readdepth+0.01 << "\t" << median_dp[i] << "\t" << l2 << endl;
		c.d_log2ratio = l2;
		lr2_values.push_back(l2);
		ofile.seekp(npos, ios::beg);
		ofile.write((char *)(&c), sizeof(c));
	}
	file.close();
	ofile.close();
return lr2_values;
}


void update_g_w_ratio (string binname, vector<double> gc_corrected, vector<double> wave_corrected, int filename_position, int number_lines) {
	rd c;
	ifstream file(binname.c_str(), ios::binary);
	ofstream ofile(binname.c_str(), ios::binary | ios::out | ios::in);
	long starting_position = filename_position*(sizeof(c)*number_lines);
	file.seekg(starting_position);
	for(int i=0; i<number_lines;i++) {
		long npos = starting_position+(sizeof(c)*i);
		file.seekg(npos);
		file.read((char *)(&c), sizeof(c));
		c.g_log2ratio = gc_corrected[i];
		c.w_log2ratio = wave_corrected[i];
		ofile.seekp(npos, ios::beg);
		ofile.write((char *)(&c), sizeof(c));
	}
	file.close();
	ofile.close();
}

vector <double> generate_sample_log2_ratio(int sample_index, string binname, vector<string> filenames, int number_lines, int max_c_filenames) {
	vector <string> csamples = get_highest_correlated_filenames(binname, number_lines, filenames, sample_index, max_c_filenames);
	vector<int> sindex = get_sample_indexes(filenames, csamples);
	vector<double> median_dp = get_median_readdepth(binname, filenames, sindex, number_lines);
	vector <double> lgr_values = update_log2_ratio(binname, filenames, median_dp, sample_index, number_lines);
return lgr_values;
}

void get_all_log2_ratios(string binname, vector<string> filenames, int filename_position, int number_lines) {
	rd c;
	ifstream file(binname.c_str(), ios::binary);
	ofstream ofile(binname.c_str(), ios::binary | ios::out | ios::in);
	long starting_position = filename_position*(sizeof(c)*number_lines);
	file.seekg(starting_position);	
	for(int j=0;j<number_lines;j++) {
		vector < rd > pp = get_probe( binname, filenames,  j, number_lines);
		vector<double> lr2_values;
		for(int i=0;i<pp.size();i++) {
			cout << pp[i].w_log2ratio << "\t";
			lr2_values.push_back(pp[i].w_log2ratio);
		}
	double local_maxima = find_maxima(lr2_values, floor(lr2_values.size()/3));
	cout << local_maxima << endl;
	long npos = starting_position+(sizeof(c)*j);
	file.seekg(npos);
	file.read((char *)(&c), sizeof(c));
	c.s_log2ratio = c.w_log2ratio-local_maxima;
	ofile.seekp(npos, ios::beg);
	ofile.write((char *)(&c), sizeof(c));
	//cout << local_maxima << endl;
	}
	file.close();
	ofile.close();
}

char* getCmdOption(char ** begin, char ** end, const std::string & option) {
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option) {
    return std::find(begin, end, option) != end;
}


int main(int argc, char *argv[]) {

	int i;
	double factor = 0.5;
	string exetype =  argv[1];
	string inputfile = argv[2];
	string bait_design = argv[3];
	char* pPath;
  	pPath = getenv ("LSB_JOBINDEX");
	i = atoi (pPath); i--;
	int sample_index = 0;
	string binname = "/lustre/scratch113/projects/ddd/users/tf2/rd_bin/bigdata.bin";
	//vector<string> filenames = get_filenames("small.txt");
	vector<string> filenames = get_filenames(inputfile);
	int number_lines = get_number_lines(filenames[0]);
	
	if(exetype == "rdbin") {
		makerdbin(filenames, binname);
	}
	
	if(exetype == "l2r") {
		int max_c_filenames = 25;
		vector<double> gc_values = getData(bait_design, 3);
		for(int i=0;i<filenames.size();i++) {
			vector<double> lg2_ratios = generate_sample_log2_ratio(i, binname, filenames, number_lines, max_c_filenames);
			vector<double> gc_corrected = gc_correct(gc_values, lg2_ratios);
			vector<double> wav = dtWaveXfm(gc_corrected, 1);
			vector<double> wave_corrected = rtab(wav,gc_corrected, factor, 1);
			update_g_w_ratio(binname, gc_corrected, wave_corrected, i, number_lines);
			cout << i << endl;
		}
	}
	if(exetype == "get_sample") {
		vector < rd > nrds =  get_sample(binname, filenames, i, number_lines);
		for(int i=0;i<nrds.size();i++) {
			cout << nrds[i].d_readcount << "\t" << nrds[i].d_readdepth << "\t" << nrds[i].d_log2ratio << "\t" << nrds[i].g_log2ratio << "\t" << nrds[i].w_log2ratio << "\t" << nrds[i].s_log2ratio << endl;
		}	
	}
	if(exetype == "get_probe") {
		vector <double> lr2_values;
		vector < rd > nrds =  get_probe(binname, filenames, i, number_lines);
		for(int i=0;i<nrds.size();i++) {
			cout << nrds[i].d_readcount << "\t" << nrds[i].d_readdepth << "\t" << nrds[i].d_log2ratio << "\t" << nrds[i].g_log2ratio << "\t" << nrds[i].w_log2ratio << endl;
			lr2_values.push_back(nrds[i].w_log2ratio);
		}	
		double local_maxima = find_maxima(lr2_values, floor(lr2_values.size()/3));
		cout << local_maxima << endl;
	}
	
	if(exetype == "cross_sample") {
		get_all_log2_ratios(binname, filenames, i, number_lines);
	}
	
/*	
	if (pPath!=NULL) {
		int i = atoi (pPath);
		int max_c_filenames = 25;
		i--;
		//printf ("The current index is: %i\n", i);
		vector<double> lg2_ratios = generate_sample_log2_ratio(i, binname, filenames, number_lines, max_c_filenames);
		vector<double> gc_corrected = gc_correct(gc_values, lg2_ratios);
		vector<double> wav = dtWaveXfm(gc_corrected, 1);
		vector<double> wave_corrected = rtab(wav,gc_corrected, 1);
		update_g_w_ratio(binname, gc_corrected, wave_corrected, i, number_lines);
		vector < rd > nrds =  get_sample(binname, filenames, i, number_lines);
		for(int i=0;i<nrds.size();i++) {
			cout << nrds[i].d_readcount << "\t" << nrds[i].d_readdepth << "\t" << nrds[i].d_log2ratio << "\t" << nrds[i].g_log2ratio << "\t" << nrds[i].w_log2ratio << endl;
		}	
	} else {
		makerdbin(filenames, binname);
	}
 */

	//int max_c_filenames = 200;
	//generate_sample_log2_ratio(i, binname, filenames, number_lines, max_c_filenames);
	
/*
	int sample_index = 0;
	string binname = "bigdata.bin";
	//vector<string> filenames = get_filenames("3.txt");
	vector<string> filenames = get_filenames("small.txt");
	int number_lines = get_number_lines(filenames[0]);

	
	makerdbin(filenames, binname);	
	int max_c_filenames = 2;
	for(int i=0;i<filenames.size();i++) {
		generate_sample_log2_ratio(i, binname, filenames, number_lines, max_c_filenames);
	}
	
	vector < rd > nrds =  get_sample(binname, filenames, sample_index, number_lines);
	for(int i=0;i<nrds.size();i++) {
		cout << nrds[i].d_readcount << "\t" << nrds[i].d_readdepth << "\t" << nrds[i].d_log2ratio << endl;
	}	
/*
	for(int j=0;j<number_lines;j++) {
		vector < rd > nrds = get_probe(binname, filenames, j, number_lines);
		for(int i=0;i<nrds.size();i++) {
			//cout << nrds[i].d_readdepth << "\t" << nrds[i].d_log2ratio << endl;
			cout << nrds[i].d_log2ratio << "\t";
		}
		cout << endl;
	}
	 */
	
return 0;
}


/*

void read_read_depth_at_position (string binname, vector<string> filenames, long starting_pos, int number_row_to_read) {

	rd c;
	ifstream file(binname.c_str(), ios::binary); 
	file.seekg(starting_pos);
	
	for(int i=0; i<filenames.size();i++) {
		file.read((char *)(&c), sizeof(c)); 
		cout << filenames[i] << "\t" << c.d_readdepth << endl;
	}
	
	file.close();

}


vector < rd > get_sample (string binname, vector<string> filenames, int filename_position, int number_lines) {
	rd c;
	vector< rd > sample;
	ifstream file(binname.c_str(), ios::binary); 
	file.seekg(filename_position);
	for(int i=0; i<number_lines;i++) {
		long npos = (filenames.size()*sizeof(c))*i;
		file.seekg(npos);
		file.read((char *)(&c), sizeof(c)); 
		cout << filenames[filename_position] << "\t" << c.d_readdepth << endl;
		sample.push_back(c);
	}
	file.close();
}

vector<string> seek_to_line(ifstream* file, long starting_pos, int number_row_to_read) {
	int count;
	string s; vector<string> items;
	(*file).seekg(starting_pos);
	for(int i=0;i<number_row_to_read;i++) {
 		std::getline(*file, s);
 		items.push_back(s);
 	}
return items;
}

vector< double > parse_double_items(vector<string> items, int index) {
	vector<double> dvalues;
	for(int i=0;i<items.size();i++) {
		string s = items[i];
		stringstream ss; ss << s;
  		vector<string> v;
  		while (getline( ss, s, '\t' )) v.push_back( s );
 		istringstream buffer(v[index]);
 		double dvalue; buffer >> dvalue;
 		dvalues.push_back(dvalue);
 	}
return dvalues;
}

// Utility function to create binary format from exome input file format 
void makerdbin(vector<string> filenames, string binname, int length) {	
	vector<long> positions;
	for(int i=0;i<filenames.size();i++) {
		positions.push_back(0);
	}
	ofstream binfile(binname.c_str(), std::ios::out | std::ios::binary);
	for(int j=0;j<length;j++) {
		for(int i=0;i<filenames.size();i++) {
			ifstream *f(new ifstream(filenames[i].c_str(), std::ios::in));
			vector<string> items = seek_to_line(f, positions[i], 1);
			positions[i] = (*f).tellg();
			(*f).close();
			vector <double> depths = parse_double_items(items, 4);
			for(int k=0;k<depths.size();k++) {
				rd rd_member(depths[k]);
				binfile.write((char *)(&rd_member), sizeof(rd_member));
			}
		}
		cout << j << endl;
	}
	binfile.close();

}

vector <ifstream *>  get_file_handles(vector<string> filenames) {
	vector <ifstream *> ifs;
  	for(int i=0;i<filenames.size();i++) {
  		ifstream *f(new ifstream(filenames[i].c_str(), std::ios::in));
		ifs.push_back(f); 
	}
return ifs;
}

vector<long> count_bytes(string filename) {
	string s; long nby=0;
	vector<long> nbytes;
	ifstream file( filename.c_str() );
	while (std::getline(file, s)) {
		nbytes.push_back(nby);
		nby+=s.length()+1;
	}
	 file.close();
return nbytes;
}

/*	
	vector< vector<long> > byte_offsets;
	for(int i=0;i<filenames1.size();i++) {
		cout << i << endl;
		cout << filenames1[i] << endl;
		vector<long> nbytes = count_bytes(ifs[i]);
		byte_offsets.push_back(nbytes);
	}
*/	

	//vector <ifstream *> ifs  = get_file_handles(filenames);
/*	long numbytes;
	vector< long > start_poses;
	vector< long > sizes;
	string binname = "tmpindex.bin";
	ofstream binfile(binname.c_str(), std::ios::out | std::ios::binary);
	for(int i=0;i<filenames.size();i++) {
		vector<long> nbytes = count_bytes(filenames[i]);
		binfile.write((char *)(&nbytes), sizeof(nbytes));
		sizes.push_back(nbytes.size());
		start_poses.push_back(numbytes);
		numbytes +=sizeof(nbytes);
		cout << i << endl;
	}
	binfile.close();
	
	vector< long > nbytes = read_number_bytes(binname, start_poses[0], sizes[0]);
*/
/*	
	cout << "Byte calculation complete" << endl;
	makerdbin(filenames1, "bigdata.bin", byte_offsets);
	cout << "Indexed Data Structure Written" << endl;
*/
//index_vec_t indices(3);
//std::generate(indices.begin(), indices.end(), SequenceGen(0));
//indices are {0, 1, 2}

/*
vector<double> indices;
for(int i=0;i<3;i++) {
	indices.push_back(i);
}

vector<double> Index; 
Index.push_back(3);
Index.push_back(1);
Index.push_back(2);

vector<string> Values;
Values.push_back("Third");
Values.push_back("First");
Values.push_back("Second");



std::sort(indices.begin(), indices.end(), Comp(Index));

for(int i=0;i<indices.size();i++) {
	cout << indices[i] << "\t" << Values[indices[i]] << endl;
}
*/

