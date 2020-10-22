#include <vector>

using namespace std;

template <typename T> 
vector<T> lm_resid(const vector<T>& xData, const vector<T>& yData){
    T xSum = 0, ySum = 0, xxSum = 0, xySum = 0, slope, intercept;
    for (long i = 0; i < xData.size(); i++) {
        xSum += xData[i];
        ySum += yData[i];
        xxSum += xData[i] * xData[i];
        xySum += xData[i] * yData[i];
    }    
    slope = (yData.size() * xySum - xSum * ySum) / (yData.size() * xxSum - xSum * xSum);
    intercept = (ySum - slope * xSum) / yData.size();
    std::vector<T> res;    
    for(int i=0;i<xData.size();i++) {
    	T predicted = intercept+slope*xData[i];
    	T residual = yData[i]-predicted;
    	res.push_back(residual);
    }
    
	return res;
}

vector<double> gc_correct (vector<double> gc_values, vector<double> ratio_values) {
	return lm_resid(gc_values, ratio_values);
}
