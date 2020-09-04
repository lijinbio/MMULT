#include <boost/algorithm/string.hpp>
#include <fstream>
#include "types.h"

using namespace std;

string double2string(double d, int precision) {
	ostringstream sstr;
	sstr << setprecision(precision) << d;
	return sstr.str();
}

string tostring(double val) {
	return double2string(val);
}

string tostring(string& val) {
	return val;
}

int getchrs(vector< vector< string > >& methfiles, vector< string >& chroms) {
	set< string > chrs;
	for (vector< string >&groupfiles : methfiles) {
		for (string& methfile : groupfiles) {
			ifstream fin(methfile);
			string line;
			while ((fin.good() && !fin.eof()) && getline(fin, line)) {
				if (line.empty() || line[0]=='#') continue;
				vector< string > fields;
				boost::split(fields, line, boost::is_any_of("\t"));
				chrs.insert(fields[0]);
			}
			fin.close();
		}
	}
	chroms.assign(chrs.begin(), chrs.end());
	return 0;
}

int readLaneToStrandSpecificHash(
		string& infile
		, map< string, map< int, cMeth > >& methPs
		, map< string, map< int, cMeth > >& methMs
		, string& chr
		)
{
	int colIdForChr = 0;
	int colIdForTotalCPs = 9;
	int colIdForMethCPs = 10;
	int colIdForTotalCMs = 12;
	int colIdForMethCMs = 13;
	int colIdForStart = 1;
	int colIdForNext = 7;

	string chrom;
	int start;
	int totalCPs;
	int methCPs;
	int totalCMs;
	int methCMs;
	char next;

	ifstream fin(infile);
	string line;
	while ((fin.good() && !fin.eof()) && getline(fin, line)) {
		if (line.empty() || line[0]=='#') continue;
		vector< string > fields;
		boost::split(fields, line, boost::is_any_of("\t"));
		chrom = fields[colIdForChr];
		if (chrom != chr) continue;
		start = stoi(fields[colIdForStart]);
		totalCPs = stoi(fields[colIdForTotalCPs]);
		methCPs = stoi(fields[colIdForMethCPs]);
		totalCMs = stoi(fields[colIdForTotalCMs]);
		methCMs = stoi(fields[colIdForMethCMs]);
		next = fields[colIdForNext][0];
		if(totalCPs > 0){ methPs[chrom][start] = cMeth(totalCPs, methCPs, '+', next);}
		if(totalCMs > 0){ methMs[chrom][start] = cMeth(totalCMs, methCMs, '-', next);}
	}
	fin.close();
	return 0;
}

int readLaneToStrandSpecificHash(
		string& infile
		, map< string, map< int, cMeth > >& methBs
		, string& chr
		)
{
	int colIdForChr = 0;
	int colIdForTotalCPs = 9;
	int colIdForMethCPs = 10;
	int colIdForTotalCMs = 12;
	int colIdForMethCMs = 13;
	int colIdForStart = 1;
	int colIdForNext = 7;

	string chrom;
	int start;
	int totalCPs;
	int methCPs;
	int totalCMs;
	int methCMs;
	char next;

	ifstream fin(infile);
	string line;
	while ((fin.good() && !fin.eof()) && getline(fin, line)) {
		if (line.empty() || line[0]=='#') continue;
		vector< string > fields;
		boost::split(fields, line, boost::is_any_of("\t"));
		chrom = fields[colIdForChr];
		if (chrom != chr) continue;
		start = stoi(fields[colIdForStart]);
		totalCPs = stoi(fields[colIdForTotalCPs]);
		methCPs = stoi(fields[colIdForMethCPs]);
		totalCMs = stoi(fields[colIdForTotalCMs]);
		methCMs = stoi(fields[colIdForMethCMs]);
		next = fields[colIdForNext][0];
		if (totalCPs+totalCMs>0) {
			methBs[chrom][start]=cMeth(totalCPs+totalCMs, methCPs+methCMs, 'B', next);
		}
	}
	fin.close();
	return 0;
}

int readLaneToStrandSpecificHash(
		string& infile
		, map< int, double >& methratios
		, string& chr
		, int mindepth
		)
{
	int colIdForChr=0;
	int colIdForStart=1;
	int colIdForRatio=3;
	int colIdForTotalC=4;

	string chrom;
	int start;
	double ratio;
	int totalC;

	ifstream fin(infile);
	string line;
	while ((fin.good() && !fin.eof()) && getline(fin, line)) {
		if (line.empty() || line[0]=='#') continue;
		vector< string > fields;
		boost::split(fields, line, boost::is_any_of("\t"));
		chrom = fields[colIdForChr];
		if (chrom != chr) continue;
		totalC=stoi(fields[colIdForTotalC]);
		if (totalC>=mindepth) {
			start=stoi(fields[colIdForStart]);
			ratio=stold(fields[colIdForRatio]);
			methratios[start]=ratio;
		}
	}
	fin.close();
	return 0;
}
