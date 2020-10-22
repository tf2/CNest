#include <vector>
#include <algorithm> 
#include <cmath>

using namespace std;

double median(vector<double> scores) {
  double median;
  size_t size = scores.size();
  sort(scores.begin(), scores.end());
  if (size  % 2 == 0) {
      median = (scores[size / 2 - 1] + scores[size / 2]) / 2;
  } else  {
      median = scores[size / 2];
  }
return median;
}

double sum (vector<double> xx) {
	double total = 0;
	for(int i=0;i<xx.size();i++) {
		total+=xx[i];
	}
return total;
}

double mad(vector<double> xx) {	
	double med = median(xx);
	for(int i = 0; i < xx.size(); i++) {
		xx[i] = abs(xx[i] - med);
	}
return 1.4826 * median(xx);
}

double mean(vector<double> xx) {
    double sum = 0;
    for(int i = 0; i < xx.size(); i++) {
        sum += xx[i];
    }
return sum / xx.size();
}

double var(vector<double> xx) {
	double m = mean(xx);
    double x = 0;
    for(int i=0 ; i<xx.size(); i++) {
    	x += (m-xx[i])*(m-xx[i]);
    }
return x/xx.size();
}

double corpearson(vector<double> x, vector<double> y) {
	double meanX = mean(x);
	double meanY = mean(y);
	double covXY = 0, varX = 0, varY = 0;
	for(int i = 0; i < x.size(); i++) {
		double tmpX = x[i] - meanX, tmpY = y[i] - meanY;  // centered values
		covXY += tmpX * tmpY;
		varX += tmpX * tmpX;
		varY += tmpY * tmpY;
	}
return covXY / (sqrt(varX) * sqrt(varY));
}

double quantile(vector<double> values, double quantile) {
    sort(values.begin(), values.end());	
return values[(int) (values.size() * quantile)];
}

double IQR(vector<double> values) {
   	double q25 = quantile(values, 0.25); 
   	double q75 = quantile(values, 0.75);
return q75 - q25;
}

double find_maxima (vector<double> data, int breaks) {
	int n = data.size();
	if(breaks > n) breaks = n;
	int space = floor (n/breaks);
	sort (data.begin(), data.end());
	double distance = ceil(fabs(data[0]-data[n-1]))/breaks; 
	double breakpoints[breaks+1];
	for(int i=0;i<breaks;i++) {
		breakpoints[i] = data[0]+(distance*i);
	}
	int mmax=0; int counts[breaks]; double mmean=0;
	breakpoints[breaks] = data[0]+(distance*(breaks));
	for(int i=0;i<breaks;i++) {
		int count = 0;
		for(int j=0;j<n;j++) {
			if(data[j]>breakpoints[i+1]) break;
			if(data[j]>breakpoints[i] & data[j]<=breakpoints[i+1]) count++;		
		}
		if(count > mmax) mmax=count;
		counts[i] = count;
	}
	for(int i=0;i<breaks;i++) {
		if(counts[i] == mmax) {
			if(mmean>0) {
				mmean = (mmean + ((breakpoints[i]+breakpoints[i+1])/2))/2;
			} else {
				mmean = (breakpoints[i]+breakpoints[i+1])/2;
			}
		}
	}
return mmean;
}

vector<double> Absol(vector<double> values) {
	vector<double> a(values.size());
	for (int i =0; i<values.size() ; i++)  {
		a[i] = fabs(values[i]); 
	}
return a;
}

vector<double> diff(vector<double> values) {
   	vector<double> diff(values.size());
    for (int i = 0; i<values.size() -1;i++) {
        int k = i+1; diff[i] = values[i] - values[k];
    }
return diff;
}

double dydLRs(vector<double> data, int space) {
	vector<double> x;
	for(int i=1;i<data.size();i+=space) {
		x.push_back(data[i]);
	}
	double w = IQR(diff(x))/(1.907745);
return w;
}

double ws (vector<double> data) {
	vector<double> vv; int limit = 200;
	for(int i=1;i<limit;i++) {
		double v = dydLRs(data, i);
		vv.push_back( v );
	}
return sum(Absol(diff(vv)));
}

vector<double> rtab(vector<double> p1, vector<double> p2, double factor, int automatic) { 
	if(automatic) {
		factor = ws(p2);
	}
	vector<double> pp(p2.size());
	double r = quantile(Absol(p1),0.68)*factor;
	for(int i=0;i<p2.size();i++) {
		if(fabs(p1[i])<r) pp[i] =  p2[i]-p1[i]; 
		else pp[i] = p2[i];
	}
return pp;
}

vector<double> extend(vector<double> data) {	
    vector<double> newdata(data.size()+2); newdata[0] = data[0];
    for (int i=0;i<data.size();i++) newdata[i+1] = data[i];
    newdata[newdata.size()-1] = data[data.size()-1];
return newdata;
}

vector<double> convolve(vector<double> data, vector<double> operat, vector<double> output){	
    int dataLen = data.size(); int operatorLen = operat.size();	
    for(int i = 0;i < dataLen-operatorLen+1;i++){
        output[i] = 0; for(int j = operatorLen-1;j >= 0;j--) output[i] += data[i+j]*operat[j];
    }
return output;
}

vector<double> Reflect(vector<double> dataX, double minx, double maxx) {	
    vector<double> dataY(dataX.size());
    for (int i=0;i<dataX.size();i++) {		
        if (dataX[i]>maxx) dataY[i] = (2* maxx) - dataX[i];
        else dataY[i] = dataX[i];
    }
    for (int i=0;i<dataY.size();i++) {		
        if (dataY[i]<minx) dataY[i] = (2* minx) - dataY[i];
        else dataY[i] = dataX[i];
        if (dataX[i]>maxx) dataY[i] = (2* maxx) - dataY[i];
        else dataY[i] = dataY[i];	
    }
return dataY;
}
