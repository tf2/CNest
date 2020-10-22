#ifndef GC_H
#define GC_H

template <typename T> 
vector<T> lm_resid(const vector<T>& xData, const vector<T>& yData);
vector<double> gc_correct(vector<double> gc_value, vector<double> ratio_values);

#endif