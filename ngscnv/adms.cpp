#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <math.h>
#include <fstream>
#include <iostream>
#include <cmath> 
#include <map>
#include <util.h>
#include <fh.h>

using namespace std;

int max_bin_size = 1000;

void correct_adm3(vector<double> values) {
	map<int, double> max_bin;
	//vector<double> values = getData(fn, col);
	vector<double> sorted_values(values.size());
	for(int i=0;i<sorted_values.size();i++) {
		sorted_values[i] = values[i];
	}
	std::sort (sorted_values.begin(), sorted_values.end()); 
	for(int i=0;i<sorted_values.size();i++) {
		int k = ceil(i/max_bin_size);
		int c = max_bin[k];
		if(values[i]>c) {
			max_bin[k] = sorted_values[i];
		}
	}
	map<double, double> umax_bin;
	for(map<int, double >::const_iterator it = max_bin.begin(); it != max_bin.end(); it++) {
		umax_bin[it->second]=0;	
	}
	vector<double> unique_bins;
	for(map<double, double >::const_iterator it = umax_bin.begin(); it != umax_bin.end(); it++) {
		unique_bins.push_back(it->first);
	}
	for(int i=0;i<values.size();i++) {
		int pin =0;double low =0; double high =0;
		while(values[i]>=unique_bins[pin] & pin <unique_bins.size()-2) {
			low = unique_bins[pin];
			pin++;
		}
		high=unique_bins[pin];
		pin = 0;
		vector <double> ratio_values;
		for(int j=0;j<sorted_values.size();j++) {
			if(sorted_values[j] >= low & sorted_values[j] <= high) {
				ratio_values.push_back(sorted_values[j]);
			}
			if(sorted_values[j]>high) 
				break;
		}
		double m = 0; 
		if(ratio_values.size()>0 ) {
			m = mad(ratio_values);
		}
		cout << i << "\t" << ratio_values.size() << "\t" << values[i] << "\t" << low << "\t" << high << "\t" << m << endl;
	}
}

/*
int main_adms(int argc, char *argv[]) {
	string fn = argv[1];
	map<int, double> max_bin;
	vector<double> values = getData(fn, 0);
	vector<double> sorted_values(values.size());
	for(int i=0;i<sorted_values.size();i++) {
		sorted_values[i] = values[i];
	}
	std::sort (sorted_values.begin(), sorted_values.end()); 
	for(int i=0;i<sorted_values.size();i++) {
		int k = ceil(i/max_bin_size);
		int c = max_bin[k];
		if(values[i]>c) {
			max_bin[k] = sorted_values[i];
		}
	}
	map<double, double> umax_bin;
	for(map<int, double >::const_iterator it = max_bin.begin(); it != max_bin.end(); it++) {
		umax_bin[it->second]=0;	
	}
	vector<double> unique_bins;
	for(map<double, double >::const_iterator it = umax_bin.begin(); it != umax_bin.end(); it++) {
		unique_bins.push_back(it->first);
	}
	for(int i=0;i<values.size();i++) {
		int pin =0;double low =0; double high =0;
		while(values[i]>=unique_bins[pin] & pin <unique_bins.size()-2) {
			low = unique_bins[pin];
			pin++;
		}
		high=unique_bins[pin];
		pin = 0;
		vector <double> ratio_values;
		for(int j=0;j<sorted_values.size();j++) {
			if(sorted_values[j] >= low & sorted_values[j] <= high) {
				ratio_values.push_back(sorted_values[j]);
			}
			if(sorted_values[j]>high) 
				break;
		}
		double m = 0; 
		if(ratio_values.size()>0 ) {
			m = mad(ratio_values);
		}
	cout << i << "\t" << ratio_values.size() << "\t" << values[i] << "\t" << low << "\t" << high << "\t" << m << endl;
	}

return 0;
}

*/