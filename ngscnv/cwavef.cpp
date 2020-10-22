#include <vector>
#include <math.h>
#include <util.h>
#include <fh.h>

using namespace std;

static double g0a[14] = {-0.0046, -0.0054, 0.0170, 0.0238, -0.1067, 0.0119, 0.5688, 0.7561, 0.2753, -0.1172, -0.0389, 0.0347, -0.0039, 0.0033};
static double g0b[14] = {0.0033, -0.0039, 0.0347, -0.0389, -0.1172, 0.2753, 0.7561, 0.5688, 0.0119, -0.1067, 0.0238, 0.0170, -0.0054, -0.0046};
static double g1a[14] = {-0.0033, -0.0039, -0.0347, -0.0389, 0.1172, 0.2753, -0.7561, 0.5688, -0.0119, -0.1067, -0.0238, 0.0170, 0.0054, -0.0046};
static double g1b[14] = {-0.0046, 0.0054, 0.0170, -0.0238, -0.1067, -0.0119, 0.5688, -0.7561, 0.2753, 0.1172, -0.0389, -0.0347, -0.0039, -0.0033};	
static double g0o[19] = {0.0001, 0, -0.0013, -0.0019, 0.0072, 0.0239, -0.0556, -0.0517, 0.2998, 0.5594, 0.2998, -0.0517, -0.0556, 0.0239, 0.0072, -0.0019, -0.0013, 0, 0.0001};
static double g1o[13] = {-0.0018, 0, 0.0223, 0.0469, -0.0482, -0.2969, 0.5555, -0.2969, -0.0482, 0.0469, 0.0223, 0, -0.0018};
static vector<double> vg0a (g0a, g0a + sizeof(g0a) / sizeof(g0a[0]) );
static vector<double> vg0b (g0b, g0b + sizeof(g0b) / sizeof(g0b[0]) );
static vector<double> vg1a (g1a, g1a + sizeof(g1a) / sizeof(g1a[0]) );
static vector<double> vg1b (g1b, g1b + sizeof(g1b) / sizeof(g1b[0]) );
static vector<double> vg0o (g0o, g0o + sizeof(g0o) / sizeof(g0o[0]) );
static vector<double> vg1o (g1o, g1o + sizeof(g1o) / sizeof(g1o[0]) );	
	
static double h0a[14] = {0.0033, -0.0039, 0.0347, -0.0389, -0.1172, 0.2753, 0.7561, 0.5688, 0.0119, -0.1067, 0.0238, 0.0170, -0.0054, -0.0046};	
static double h0b[14] = {-0.0046, -0.0054, 0.0170, 0.0238, -0.1067, 0.0119, 0.5688, 0.7561, 0.2753, -0.1172, -0.0389, 0.0347, -0.0039, 0.0033};	
static double h1a[14] = {-0.0046, 0.0054, 0.0170, -0.0238, -0.1067, -0.0119, 0.5688, -0.7561, 0.2753, 0.1172, -0.0389, -0.0347, -0.0039, -0.0033};
static double h1b[14] = {-0.0033, -0.0039, -0.0347, -0.0389, 0.1172, 0.2753, -0.7561, 0.5688, -0.0119, -0.1067, -0.0238, 0.0170, 0.0054, -0.0046};
static double h0o[13] = {-0.0018, 0, 0.0223, -0.0469, -0.0482, 0.2969, 0.5555, 0.2969, -0.0482, -0.0469, 0.0223, 0, -0.0018};
static double h1o[19] = {-0.0001, 0, 0.0013, -0.0019, -0.0072, 0.0239, 0.0556, -0.0517, -0.2998, 0.5594, -0.2998, -0.0517, 0.0556, 0.0239, -0.0072, -0.0019, 0.0013, 0, -0.0001};
static vector<double> vh0a (h0a, h0a + sizeof(h0a) / sizeof(h0a[0]) );
static vector<double> vh0b (h0b, h0b + sizeof(h0b) / sizeof(h0b[0]) );
static vector<double> vh1a (h1a, h1a + sizeof(h1a) / sizeof(h1a[0]) );
static vector<double> vh1b (h1b, h1b + sizeof(h1b) / sizeof(h1b[0]) );
static vector<double> vh0o (h0o, h0o + sizeof(h0o) / sizeof(h0o[0]) );
static vector<double> vh1o (h1o, h1o + sizeof(h1o) / sizeof(h1o[0]) );


vector<double> colFilter(vector<double> data, vector<double> filter) {	
    int r = data.size(); 
    int m = filter.size(); 
    int m2 = (int)floor(m/2);	
    vector<double> reF;
    for (int i=1-m2; i<r+m2;i++) { reF.push_back(i); }  
    vector<double> res = Reflect(reF, 0.05, r+0.05);
	for (int i=0;i<res.size();i++) res[i] = (int)ceil(res[i]);
    vector<double> res2(res.size());
    for (int i=0;i<res2.size();i++) { int pin2 = (int) res[i]-1; res2[i] = data[pin2]; }	
    vector<double> b(res2.size());
    vector<double> c = convolve(res2, filter, b);
return c;	
}

vector<double> coldFilt(vector<double> X, vector<double> ha, vector<double> hb) {		
    int r = X.size(); int m = ha.size(); int m2 = (int)floor(m/2);		
    vector<double> reF;
    for (int i=1-m; i<r+m+1;i++) { reF.push_back(i); }
    vector<double> reF2;
    for (int i=0;i<reF.size();i++) reF2.push_back(reF[i]);
    vector<double> res = Reflect(reF2, 0.05, r+0.05);	
    for (int i=0;i<res.size();i++) res[i] = (int)ceil(res[i]);
    vector<double> hao(m2); vector<double> hae(m2); vector<double> hbo(m2); vector<double> hbe(m2);		
    int pin1 =0; int pin2=1;		
    for (int i=0;i<m2;i++) {
        hao[i] = ha[pin1]; hae[i] = ha[pin2]; hbo[i] = hb[pin1]; hbe[i] = hb[pin2];
        pin1+=2; pin2+=2;
    }    
    vector<int> t;
    int pin3 = 6;
    for (int i=6;pin3<r+2*m-2;i++) { t.push_back(pin3); pin3 = pin3 +4; }
    t.push_back(pin3);
	int r2 = r/2; vector<double> Y(r2);	
    double sum =0;	
    for (int i=0;i<ha.size();i++)  sum+=ha[i]*hb[i];		
    vector<double> s1; vector<double> s2 ;
    pin1 = 1;	
    if (sum>0) {	
        for (int i=0;i<r2/2;i++) { s1.push_back(pin1); s2.push_back(pin1+1); pin1+=2; }	
    } else {
        for (int i=0;i<r2/2;i++) { s2.push_back(pin1); s1.push_back(pin1+1); pin1+=2; }
    }		
    vector<double> Fin1; vector<double> Fin2; vector<double> Fin3; vector<double> Fin4;
    for (int i=0;i<t.size();i++) {			
        int pp = t[i];    
        int ppp = (int) res[pp-2]; if (ppp==X.size()) { ppp = X.size()-1; } if (ppp<0) { ppp=0; }
        Fin1.push_back(X[ppp-1]);
        ppp = (int) res[pp-4]; if (ppp==X.size()) { ppp = X.size()-1; } if (ppp<0) { ppp=0; }
        Fin2.push_back(X[ppp-1]);
        ppp = (int) res[pp-1]; if (ppp==X.size()) { ppp = X.size()-1; } if (ppp<0) { ppp=0; }
        Fin3.push_back(X[ppp-1]);
        ppp = (int) res[pp-3]; if (ppp==X.size()) { ppp = X.size()-1; } if (ppp<0) { ppp=0; }
        Fin4.push_back(X[ppp-1]);	
    }
    vector<double> blank1(Fin1.size()); vector<double> blank2(Fin2.size()); vector<double> blank3(Fin3.size()); vector<double> blank4(Fin4.size());
    vector<double> temp1 =  convolve(Fin1,hao,blank1); 
    vector<double> temp2 =  convolve(Fin2,hae,blank2);
    vector<double> temp3 =  convolve(Fin3,hbo,blank3); 
    vector<double> temp4 =  convolve(Fin4,hbe,blank4);				
    for (int i=0;i<r2/2;i++) {
        int p1 = (int) s2[i]-1; int p2 = (int) s1[i]-1;
        Y[p1] = temp1[i] + temp2[i]; Y[p2] = temp3[i] + temp4[i];       
    }	
return Y;		
}

	
vector<double> coliFilt(vector<double> X, vector<double> ha, vector<double> hb) {
        int r = X.size(); int m = ha.size(); int m2 = (int) floor(m/2); vector<double> Y(r*2);
        vector<double> hao(m2); vector<double> hae(m2);
        vector<double> hbo(m2); vector<double> hbe(m2);	
        int pin1 =0; int pin2=1;		
        for (int i=0;i<m2;i++) {
            hao[i] = ha[pin1]; hae[i] = ha[pin2];		
            hbo[i] = hb[pin1]; hbe[i] = hb[pin2];
            pin1+=2; pin2+=2;
        }
        vector<double> reF;	
        int pin =0;
        for (int i=1-m2; i<r+m2+1;i++) { reF.push_back(i); pin++;  }		
        vector<double> reF2(reF.size());	
        for (int i=0;i<reF.size();i++) reF2[i] = reF[i];	
        vector<double> res = Reflect(reF2, 0.05, r+0.05);
        for (int i=0;i<res.size();i++) res[i] = (int)ceil(res[i]);		
        vector<double> res2(res.size());	
        for (int i=0;i<res.size();i++) { int pin3 = (int) res[i]-1; res2[i] = X[pin3]; }	
        if (m2 % 2 == 0) {		
			vector<double> t; int ppin = 4;
			for (int i=4;i<m+r/2-4;i++) { t.push_back(ppin); ppin=ppin+2; }		
			double sum =0;
			for (int i=0;i<ha.size();i++) sum+=ha[i]*hb[i];
			vector<double> ta(t.size()); vector<double> tb(t.size());
			
			if (sum>0) {	
				for (int i=0;i<t.size();i++) { ta[i] = t[i]; tb[i] = -1; }			
			} else {
                for (int i=0;i<t.size();i++) { ta[i] = t[i]-1; tb[i] = t[i]; }
			}		
			vector<double> t1(ta.size()); vector<double> t2(ta.size());			
			vector<double> t3(ta.size()); vector<double> t4(ta.size());			
			for (int i=0;i<ta.size();i++) {	
				int p1 = (int) ta[i]; int p2 = (int) tb[i];		
				t1[i] = res2[p1-1]; t2[i] = res2[p2-1];	
				int p3 = (int) ta[i]-2; int p4 = (int) tb[i]-2;	
				t3[i] = res2[p2-1]; t4[i] = res2[p3-1];	
			}
			vector<double> blank1(t4.size()); vector<double> blank2(t3.size()); vector<double> blank3(t2.size()); vector<double> blank4(t1.size());
			vector<double> temp1 =  convolve(t4,hae,blank1);
			vector<double> temp2 =  convolve(t3,hbe,blank2);
			vector<double> temp3 =  convolve(t2,hao,blank3);
			vector<double> temp4 =  convolve(t1,hbo,blank4);	
			int mPin = 0;
			for (int i=0;i<Y.size()/4;i++) {
				Y[mPin] = temp1[i]; Y[mPin+1] = temp2[i]; Y[mPin+2] = temp3[i]; Y[mPin+3] = temp4[i]; 
				mPin = mPin+4;
			}	
			
        } else {			
			vector<double> t; int ppin = 3;
			for (int i=4;i<m+r/2-5;i++) { t.push_back(ppin); ppin=ppin+2; }		
			double sum =0;
			for (int i=0;i<ha.size();i++)  sum+=ha[i]*hb[i];
			vector<double> ta(t.size()); vector<double> tb(t.size());		
			if (sum>0) {		
				for (int i=0;i<t.size();i++) { ta[i] = t[i]; tb[i] = t[i]-1; }		
			} else {
                for (int i=0;i<t.size();i++) { ta[i] = t[i]-1; tb[i] = t[i]; }
			}			
			vector<double> t1(ta.size()); vector<double> t2(ta.size());	
			vector<double> t3(ta.size()); vector<double> t4(ta.size());	
			for (int i=0;i<ta.size();i++) {				
				int p1 = (int) ta[i]; int p2 = (int) tb[i];			
				t1[i] = res2[p1-1]; t2[i] = res2[p2-1];	
				int p3 = (int) ta[i]-2; int p4 = (int) tb[i]-2;		
				t3[i] = res2[p2-1]; t4[i] = res2[p3-1];
			}		
			vector<double> blank12(t2.size()); vector<double> blank22(t1.size()); vector<double> blank32(t2.size()); vector<double> blank42(t1.size());
			vector<double> temp1 =  convolve(t2,hao,blank12);
			vector<double> temp2 =  convolve(t1,hbo,blank22);
			vector<double> temp3 =  convolve(t2,hae,blank32);
			vector<double> temp4 =  convolve(t1,hbe,blank42);			
			int mPin = 0;
			for (int i=0;i<Y.size()/4;i++) {
				Y[mPin] = temp1[i]; Y[mPin+1] = temp2[i]; Y[mPin+2] = temp3[i]; Y[mPin+3] = temp4[i];
				mPin = mPin+4;
			}
        }	
return Y;	
}

vector<double> dtWaveIfm(vector<double> Yl, vector<double> Yh, int nLevels) {		
	vector<double> Lo(Yl.size());
	for(int i=0;i<Yl.size();i++) {
		Lo[i] = Yl[i];
	}
	for (int i =0;i<nLevels;i++) Lo = coliFilt(Yl, vg0b, vg0a);	
return Lo;
}
	
vector<double> dtWaveXfm(vector<double> data, int nLevels) {		
	int L = data.size(); bool odd = false;
    if (L % 2 == 1 ) { odd=true; L++; }	
	vector<double> data2;
	for(int i=0;i<L-1;i++) data2.push_back(data[i]);
	if(odd) data2[L-1] = 0;		
    vector<double> Hi = colFilter(data2, vh1o); 
    vector<double> Lo = colFilter(data2, vh0o);  
    for (int i=0;i<nLevels;i++) {	
        if (Lo.size()<10) break;	
        if (Lo.size() % 4 != 0 ) Lo = extend(Lo);
        Hi = coldFilt(Lo,vh1b,vh1a); 
        Lo = coldFilt(Lo,vh0b,vh0a);
		Lo = dtWaveIfm(Lo, Hi, nLevels);	
    }      
    for(int i=Lo.size()-1;i<0;i++) Lo[i]=Lo[i-1];
return Lo;	
}

int main_cwavef(int argc, char *argv[]) {
	double factor = 0.5;
	int find_level = 0;
	string fn = argv[1];
    if(argc > 3) {
    	find_level = atoi(argv[2]);
    }
	vector<double> ddd = getData(fn, 3);
	vector<double> wav = dtWaveXfm(ddd, 1);
	vector<double> ws1 = rtab(wav,ddd, factor, find_level);
	for (int i=0;i<ddd.size();i++) {
		//cout << ddd[i] << "\t" << wav[i] << "\t" << ws1[i] << endl;
	}
return 0;
}