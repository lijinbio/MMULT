#ifndef __TYPES__
#define __TYPES__

#define BATCH 20000

#include <string>
#include <vector>
#include <set>
#include <map>
#include <iomanip>
#include <sstream>
#include <iostream>

using namespace std;

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
		vector< vector< string > >& methfiles
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
		);

int readLaneToStrandSpecificHash(
		string& infile
		, map< int, double >& methratios
		, string& chr
		, int mindepth=3
		);

class Region {
	public:
		int start;
		int end;
		int b; // number of CpG
		int N; // number of samples
		vector< int > cpgs;
		vector< vector< double > > ratios; // N by b matrix
	public:
		double nominal; // norminal mean methylation ratio
		double sd; // average standard deviation of methylation ratio
		double NME; // normalized methylation entropy
	public:
		void out() {
			cout << "Region: " << start << "\t" << end << "\t" << N << "\t" << b << "\t" << strjoin(cpgs) << endl;
			for (int r=0; r<N; r++) {
				cout << "row" << r+1 << ": " << strjoin(ratios[r]) << endl;
			}
		}
};

class VMC {
	public:
		string chrom;
		int start;
		int end;
		vector< double > ratios;
		int N;
	public:
		double nominal;
		double sd;
		double NME;
	public:
		void out() {
			cout << "CpG: "
				<< chrom << "\t"
				<< start << "\t"
				<< end << "\t"
				<< N << "\t"
				<< nominal << "\t"
				<< sd << "\t"
				<< NME << "\t"
				<< strjoin(ratios) << endl;
		}
};

class VMR {
	public:
		string chrom;
		int start;
		int end;
		vector< VMC > cpgs;
		int N;
	public:
		double nominal;
		double sd;
		double NME;
	public:
		void cpgs2nominal() {
			double sum=0;
			for (VMC& vmc: cpgs) {
				sum+=vmc.nominal;
			}
			nominal=sum*1.0/cpgs.size();
		}
		void cpgs2sd() {
			double sum=0;
			for (VMC& vmc: cpgs) {
				sum+=vmc.sd;
			}
			sd=sum*1.0/cpgs.size();
		}
		void cpgs2nme() {
			double sum=0;
			for (VMC& vmc: cpgs) {
				sum+=vmc.NME;
			}
			NME=sum*1.0/cpgs.size();
		}
	public:
		void out() {
			cout << "VMR: "
				<< chrom << "\t"
				<< start << "\t"
				<< end << "\t"
				<< N << "\t"
				<< nominal << "\t"
				<< sd << "\t"
				<< NME << endl;
			for (VMC& vmc: cpgs) {
				vmc.out();
			}
		}
};

const string discretizelabel="0123456789ABCDEFGHIJKLMNOPTRSTUVWXYZabcdefghijklmnoptrstuvwxyz";

#endif
