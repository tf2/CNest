#ifndef FH_H
#define FH_H

#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

using namespace std;

int get_number_lines (string filename);
vector< string > get_filenames (string filename);
vector<double> getData(string filename, int index);
vector< vector<double> > getData(string filename, vector<int> col_indexes);

#endif