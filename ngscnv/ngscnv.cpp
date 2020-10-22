#include <vector>
#include <iostream>
#include <util.h>
#include <fh.h>
#include <cwavef.h>
#include <gc.h>
#include <adms.h>
#include <readdata.h>

using namespace std;

int main(int argc, char *argv[]) {

	string project_name = std::string(argv[1]);
	string execution_type = std::string(argv[2]);
	
	if(execution_type == "new-project") {
		string command = "mkdir " + project_name;
		system(command.c_str());
	}
	
	if(execution_type == "add-baits") {
		string bait_design = std::string(argv[3]);
		string command = "cp " + bait_design + " " + project_name + "/index.txt" ;
		system(command.c_str());
	}
	
	if(execution_type == "bam-to-rd") {
		string bam_file = std::string(argv[3]);
		string sample_id = std::string(argv[4]);
		string baitfile = project_name + "/index.txt";
		string binary_outfile = project_name + "/" + sample_id;
		write_bam_to_rd(strdup(bam_file.c_str()), strdup(baitfile.c_str()), strdup(binary_outfile.c_str()));
	}
	
	if(execution_type == "rd-read") {
		string sample_id = std::string(argv[3]);
		string chr = std::string(argv[4]);
		int start = atoi(argv[5]);
		int stop = atoi(argv[6]);
		string baitfile = project_name + "/index.txt";
		string binary_outfile = project_name + "/" + sample_id;
		
		char *indexfile = strdup(baitfile.c_str());
		char *binname = strdup(binary_outfile.c_str());
		
		//vector < int >  offsets_and_rows = getOffsets(indexfile, "1", 20503, 31839);
	
		vector < vector<int> > positions = findOffset(indexfile, chr, start, stop);
		vector<int> v = positions[positions.size()-1];
		int offset = v[0]; 
		int number_rows = v[1];
		vector < readdata > data =  getValues(binname, offset, number_rows);
		for(int i=0;i<data.size();i++) {
			//cout << data[i].d_readcount << "\t" << data[i].d_depth << endl;
			vector <int> coordinates = positions[i];
			cout << coordinates[0] << ":" << coordinates[1] << "-" << coordinates[2] << "\t" << data[i].d_readcount << "\t" << data[i].d_depth << endl;
		}
	}

	if(execution_type == "rd-dump") {
		string sample_id = std::string(argv[3]);
		string binary_outfile = project_name + "/" + sample_id;
		vector<int> col_indexes(2); col_indexes[0] = 0; col_indexes[1]= 1;
        string baitfile = project_name + "/index.txt";
		char *binname = strdup(binary_outfile.c_str());
        //vector< vector<double> > values = getData(baitfile, col_indexes);
        int number_lines = get_number_lines(baitfile);
		//vector < readdata > data =  getValues(binname, 0, values[0].size());
        vector < readdata > data =  getValues(binname, 0, number_lines);
        for(int i=0;i<data.size();i++) {
            cout << data[i].d_readcount << "\t" << data[i].d_depth << endl;
		}
	
	}
	
	if(execution_type=="normalisation") {
		int find_wave_level = 0;
		double factor = 0.7;
		string fn = argv[1];
		vector<int> col_indexes(2); col_indexes[0] = 0; col_indexes[1]= 1;
		vector< vector<double> > values = getData(fn, col_indexes);
		vector<double> residuals = lm_resid(values[0], values[1]);
		vector<double> wav = dtWaveXfm(residuals, 1);
		vector<double> ws1 = rtab(wav,residuals, factor, find_wave_level);
		//correct_adm3(values[1]);
		for (int i=0;i<residuals.size();i++) {
			cout << values[0][i] << "\t" << values[1][i] << "\t" << residuals[i] << "\t" << wav[i] << "\t" << ws1[i] << endl;
		}
	}
 
return 0;
}
