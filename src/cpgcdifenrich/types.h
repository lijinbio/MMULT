#ifndef __TYPES__
#define __TYPES__

#ifdef __APPLE__
#define BATCH 2
#define UNITSIZE 10240 // Limited shared memory on Mac OS X
#else
#define BATCH 2000
#define UNITSIZE 102400
#endif

#include <string>
#include <vector>
#include <set>
#include <map>
#include <iomanip>
#include <sstream>
#include <iostream>

using namespace std;

double string2double(const string& str, int precision=3);

string double2string(double d, int precision=3);
template<class T>
string tostring(T val) {
	return std::to_string(val);
}
string tostring(double val);
string tostring(string& val);
template<class T>
string strjoin(vector< T >& v, string separator=",") {
	string s=tostring(v[0]);
	for (int i=1; i<v.size(); i++) {
		s+=separator;
		s+=tostring(v[i]);
	}
	return s;
}

template<class T, class M>
string strjoin(vector< T >& left, vector< M >& right, string fieldseparator="", string separator=",") {
	string s=tostring(left[0])+fieldseparator+tostring(right[0]);
	for (int i=1; i<left.size(); i++) {
		s+=separator;
		s+=tostring(left[i])+fieldseparator+tostring(right[i]);
	}
	return s;
}

// mcomp.cpp in MOABS
class cMeth
{
	public:
		int totalC;
		int methC;
		char strand;
		char nextN;
		cMeth(): totalC(0), methC(0), strand('*'), nextN('N') { }
		cMeth(int t, int m, char s, char n) : totalC(t), methC(m), strand(s), nextN(n) {}
};

int getchrs(
		vector< vector< string > >& infiles
		, vector< string >& chroms
		);
int getchrs(
		vector< string >& infiles
		, vector< string >& chroms
		);
// mcomp.cpp in MOABS
int readLaneToStrandSpecificHash(
		string& infile
		, map< string, map< int, cMeth > >& methPs
		, map< string, map< int, cMeth > >& methMs
		, string& chr
		);

int readLaneToStrandSpecificHash(
		string& infile
		, map< string, map< int, cMeth > >& methBs
		, string& chr
		, int mindepth=1
		);

int readLaneToStrandSpecificHash(
		string& infile
		, map< string, map< int, cMeth > >& methBs
		, string& chr
		, int lower
		, int upper
		, int mindepth=1
		);

int readLaneToStrandSpecificHash(
		string& infile
		, map< int, double >& methratios
		, string& chr
		, int mindepth=1
		);

int readCDIF(
		string& infile
		, map< int, double >& cdifs
		, string& chr
		);
int readCDIF(
		string& infile
		, map< int, double >& cdifs
		, string& chr
		, int lower
		, int upper
		);

struct BBGLMLRT {
	int sJx;
	int sKx;
	vector< int > groupid;
	vector< int > groupsize;
	vector< int > tcs;
	vector< int > mcs;
	double ssx; // MLE results
	double pvalue;
	vector< double > pi;
};

class Hist {
	public:
		vector< long long > counts;
		double lower;
		double upper;
		int numbins;
	public:
		vector< double > density;
		vector< double > pi;
		vector< string > labels;
	public:
		int counts2hist();
		int output(ostream& outfile);
		int output(string& outfile);
};

class DMC {
	public:
		string chrom;
		int start;
		int end;
		double meancdif;
		double sumcdif;
		double sd;
		double kld;
		double entropy;
		string label;
};

class DMR {
	public:
		string chrom;
		int start;
		int end;
		vector< double > meancdif;
		vector< double > kld;
		vector< double > sd;
		string label;
	public:
		DMR(): start(-1), end(-1) {}
		int clear();
};

bool pairgreater(pair< double, int >& a, pair< double, int >& b);
int padjustfdr(vector< double >& p, vector< double >& q);

#endif
