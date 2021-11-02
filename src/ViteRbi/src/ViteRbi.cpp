#include<set>
#include<map>
#include<math.h>
#include<vector>
#include<sstream>
#include<fstream>
#include<iostream>

using namespace std;

#include<R.h>
#include<Rdefines.h>
#include<Rinternals.h>
#include<Rinterface.h>

/*
Description:	HMM Intepretation.
Author:			Tomas William Fitzgerald
Date:			30/01/2012
Licence:		Artistic
*/

extern "C" {

const double Pi = 3.141592;
map < string, vector<string> > mdata;

class hmm {
	public:
    	hmm(): d_nsts(0), d_nem(0), d_ntr(0), d_emsp(0), d_trsp(0) {}
    	hmm(double nsts, int nem, int ntr, double *emsp, double *trsp): d_nsts(nsts), d_nem(nem), d_ntr(ntr), d_emsp(emsp), d_trsp(trsp){}
	public:
   		int d_nsts; int d_nem; int d_ntr; double *d_emsp; double *d_trsp;
};

double npdf(double  mu, double sigma, double x) {
	double f = (x - mu)/sigma; double p = exp(-0.5*pow(f,2) ) / (sqrt(2*Pi) * sigma);
return p;
}

vector < double** > ets(hmm h) {
	double **ep; double **tp; 
	int N = (int)h.d_nsts; ep = new double*[N]; tp = new double*[N]; int c=0;
	for (int i=0;i<N;++i) { ep[i] = new double[2]; tp[i] = new double[3]; }
	for(int i=0;i<N;i++) {
		for(int j=0;j<h.d_nem;j++) { ep[i][j]=h.d_emsp[c]; c++; }
	} c=0;
	for(int i=0;i<N;i++) {
		for(int j=0;j<h.d_ntr;j++) { tp[i][j]=h.d_trsp[c]; c++; }
	}
	vector<  double** > r; r.push_back(ep); r.push_back(tp);
return r;
}

vector<double> viterbi(hmm h, vector<double> d, vector<double> s) {        
	int i, j; int N = (int)h.d_nsts; 
	int **fwds = new int*[N]; int **bwds = new int*[N]; 
	vector < double** > rrr = ets(h); double **ep = rrr[0]; double **tp = rrr[1];
    for(i=0;i<N;++i) { fwds[i] = new int[d.size()]; bwds[i] = new int[d.size()]; }
    for(i=0;i<N;i++) { fwds[i][0] = log(npdf(ep[i][0], ep[i][1], d[0])); }  
        for (int i=1;i<d.size();i++) {
            for (int j=0;j<N;j++) {
                int pps; double mp; double sps[N];
                for (int ps = 0; ps < N; ps++) { sps[ps] = fwds[ps][i-1] + log(tp[ps][j]); }
                mp = sps[0];  pps = 0;
				for (int k=0;k<N;k++) {
					if (sps[k] > mp) { mp = sps[k]; pps = k; }
				}
                fwds[j][i] = mp + log(npdf(ep[j][0], ep[j][1], d[i])); bwds[j][i] = pps;
         	}
        }  
        int fs=0; double ml=fwds[0][d.size()-1];
        if(fwds[1][d.size()-1]>ml){ ml=fwds[1][d.size()-1]; fs=1; }
        if(fwds[2][d.size()-1]>ml){ ml=fwds[2][d.size()-1]; fs=2; }
        s[d.size()-1]=fs; int ns=fs;
        for(int i= d.size()-2;i>=0;i--){ s[i]=bwds[ns][i+1]; ns=s[i]; }
return s;        
}

vector<double> nord(double *data, int N) {
	double su=0; int count=0; vector < double > e;
	for(int i=0;i<N;i++) { su+=data[i]; count++; }
	double sd; double sus=0; double m=su/(count);    
    for (int i=0;i<count; i++) {  sus+=pow((data[i]-m),2); }
    double var = sus/(count-1); sd=sqrt(var);
    for (int i=0; i<count;i++) { data[i]=(data[i]- m)/(float)sd; e.push_back(data[i]); }
return e;
}

void ViteRbi(double *data, double *states, double *emissions, double *transitions, int *dN, int *sN, int *eN, int *tN) {
	vector < double > d = nord(data, (int)*dN);
	vector < double > s; for(int i=0;i<(int)*dN;i++) { s.push_back(states[i]); }
	hmm h((int)*sN, (int)*eN, (int)*tN, emissions, transitions);
	vector< double > r = viterbi(h,d,s);
	for(int i=0;i<(int)*dN;i++) {
		states[i]=r[i];
	}
}


}
