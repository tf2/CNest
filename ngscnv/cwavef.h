#ifndef CWAVEF_H
#define CWAVEF_H

vector<double> colFilter(vector<double> data, vector<double> filter) ;
vector<double> coldFilt(vector<double> X, vector<double> ha, vector<double> hb) ;
vector<double> coliFilt(vector<double> X, vector<double> ha, vector<double> hb) ;
vector<double> dtWaveIfm(vector<double> Yl, vector<double> Yh, int nLevels);
vector<double> dtWaveXfm(vector<double> data, int nLevels) ;
double dydLRs(vector<double> data, int space);
double ws (vector<double> data);

#endif
