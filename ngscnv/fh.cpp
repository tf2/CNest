#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

using namespace std;

int get_number_lines (string filename) {
	int number_of_lines = 0;
    std::string line;
    std::ifstream myfile(filename.c_str());
    while (std::getline(myfile, line))
        ++number_of_lines;
return number_of_lines;
}

vector< string > get_filenames (string filename) {
	ifstream file( filename.c_str() );
	string s; vector< string > filenames;
	while (std::getline(file, s)) {
 		filenames.push_back(s);
	}
	file.close();
	return filenames;
}

vector<double> getData(string filename, int index) {	
	string s; vector<double> ratio_values;
	std::ifstream file(filename.c_str());
	while (std::getline(file, s)) {
		stringstream ss; ss << s; vector<string> v;
  		while (getline( ss, s, '\t' )) v.push_back( s );
 		istringstream buffer1(v[index]); double rvalue; buffer1 >> rvalue;
 		ratio_values.push_back(rvalue);
	}
	file.close();
return ratio_values;
}

vector< vector<double> > getData(string filename, vector<int> col_indexes) {
    string s; vector< vector<double> > ratio_values;
	vector<double> r;
	for(int i=0;i<col_indexes.size();i++) {
		ratio_values.push_back(r);
        std::cout << i << std::endl;
	}
    std::cout << filename.c_str() << std::endl;
	ifstream file(filename.c_str());
	while (getline(file, s)) {
        std::cout << s << std::endl;
		stringstream ss; ss << s; vector<string> v;
  		while (getline( ss, s, '\t' )) v.push_back( s );
 		for(int i=0;i<col_indexes.size();i++) {
			istringstream buffer1(v[col_indexes[i]]); double rvalue; buffer1 >> rvalue;
			ratio_values[i].push_back(rvalue);
		}
	}
	file.close();
return ratio_values;
}
