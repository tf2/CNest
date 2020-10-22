#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <algorithm> 

using namespace std;

double median(vector<double> scores); 
double sum (vector<double> xx); 
double mad(vector<double> xx); 
double mean(vector<double> xx); 
double var(vector<double> xx);
double corpearson(vector<double> x, vector<double> y);
double quantile(vector<double> values, double quantile); 
double IQR(vector<double> values); 
double dLRs(vector<double> values);
double dydLRs(vector<double> data, int space);
double ws(vector<double> data);
double find_maxima (vector<double> data, int breaks);

vector<double> Absol(vector<double> values);
vector<double> diff(vector<double> values);
vector<double> rtab(vector<double> p1, vector<double> p2, double factor, int automatic);
vector<double> extend(vector<double> data);
vector<double> convolve(vector<double> data, vector<double> operat, vector<double> output);
vector<double> Reflect(vector<double> dataX, double minx, double maxx);

#endif
