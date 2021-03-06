#include <boost/algorithm/string.hpp>
#include <fstream>
#include <utility>
#include "types.h"

using namespace std;

double string2double(const string& str, int precision)
{
	stringstream sstr;
	sstr << setprecision(precision) << fixed << str << endl;
	double value;
	sstr >> value;
	return value;
}

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

int getchrs(vector< vector< string > >& infiles, vector< string >& chroms) {
	set< string > chrs;
	for (vector< string >&groupfiles : infiles) {
		for (string& infile : groupfiles) {
			ifstream fin(infile);
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

int getchrs(vector< string >& infiles, vector< string >& chroms) {
	set< string > chrs;
	for (string& infile: infiles) {
		ifstream fin(infile);
		string line;
		while ((fin.good() && !fin.eof()) && getline(fin, line)) {
			if (line.empty() || line[0]=='#') continue;
			vector< string > fields;
			boost::split(fields, line, boost::is_any_of("\t"));
			chrs.insert(fields[0]);
		}
		fin.close();
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
		else { // found the chr
			do {
				if (line.empty() || line[0]=='#') continue;
				vector< string > chrfields;
				boost::split(chrfields, line, boost::is_any_of("\t"));
				chrom = chrfields[colIdForChr];
				if (chrom != chr) break;
				start = stoi(chrfields[colIdForStart]);
				totalCPs = stoi(chrfields[colIdForTotalCPs]);
				methCPs = stoi(chrfields[colIdForMethCPs]);
				totalCMs = stoi(chrfields[colIdForTotalCMs]);
				methCMs = stoi(chrfields[colIdForMethCMs]);
				next = chrfields[colIdForNext][0];
				if(totalCPs > 0){ methPs[chrom][start] = cMeth(totalCPs, methCPs, '+', next);}
				if(totalCMs > 0){ methMs[chrom][start] = cMeth(totalCMs, methCMs, '-', next);}
			} while ((fin.good() && !fin.eof()) && getline(fin, line));
			break;
		}
	}
	fin.close();
	return 0;
}

int readLaneToStrandSpecificHash(
		string& infile
		, map< string, map< int, cMeth > >& methBs
		, string& chr
		, int mindepth
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
		else { // found the chr
			do {
				if (line.empty() || line[0]=='#') continue;
				vector< string > chrfields;
				boost::split(chrfields, line, boost::is_any_of("\t"));
				chrom = chrfields[colIdForChr];
				if (chrom != chr) break;
				start = stoi(chrfields[colIdForStart]);
				totalCPs = stoi(chrfields[colIdForTotalCPs]);
				methCPs = stoi(chrfields[colIdForMethCPs]);
				totalCMs = stoi(chrfields[colIdForTotalCMs]);
				methCMs = stoi(chrfields[colIdForMethCMs]);
				next = chrfields[colIdForNext][0];
				if (totalCPs+totalCMs>=mindepth) {
					methBs[chrom][start]=cMeth(totalCPs+totalCMs, methCPs+methCMs, 'B', next);
				}
			} while ((fin.good() && !fin.eof()) && getline(fin, line));
			break;
		}
	}
	fin.close();
	return 0;
}

int readLaneToStrandSpecificHash(
		string& infile
		, map< string, map< int, cMeth > >& methBs
		, string& chr
		, int lower
		, int upper
		, int mindepth
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
		else { // found the chr
			do {
				if (line.empty() || line[0]=='#') continue;
				vector< string > chrfields;
				boost::split(chrfields, line, boost::is_any_of("\t"));
				chrom = chrfields[colIdForChr];
				if (chrom != chr) break;
				start = stoi(chrfields[colIdForStart]);
				if (start<lower) continue; // [lower, upper)
				if (start>=upper) break; // BED file should be sorted within a chromosome.
				totalCPs = stoi(chrfields[colIdForTotalCPs]);
				methCPs = stoi(chrfields[colIdForMethCPs]);
				totalCMs = stoi(chrfields[colIdForTotalCMs]);
				methCMs = stoi(chrfields[colIdForMethCMs]);
				next = chrfields[colIdForNext][0];
				if (totalCPs+totalCMs>=mindepth) {
					methBs[chrom][start]=cMeth(totalCPs+totalCMs, methCPs+methCMs, 'B', next);
				}
			} while ((fin.good() && !fin.eof()) && getline(fin, line));
			break;
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
		else {
			do {
				if (line.empty() || line[0]=='#') continue;
				vector< string > chrfields;
				boost::split(chrfields, line, boost::is_any_of("\t"));
				chrom = chrfields[colIdForChr];
				if (chrom != chr) break;
				totalC=stoi(chrfields[colIdForTotalC]);
				if (totalC>=mindepth) {
					start=stoi(chrfields[colIdForStart]);
					ratio=stold(chrfields[colIdForRatio]);
					methratios[start]=ratio;
				}
			} while ((fin.good() && !fin.eof()) && getline(fin, line));
			break;
		}
	}
	fin.close();
	return 0;
}

int readCDIF(
		string& infile
		, map< int, double >& cdifs
		, string& chr
		)
{
	int colIdForChr=0;
	int colIdForStart=1;
	int colIdForCDIF=13;

	string chrom;
	int start;
	double cdif;

	ifstream fin(infile);
	string line;
	while ((fin.good() && !fin.eof()) && getline(fin, line)) {
		if (line.empty() || line[0]=='#') continue;
		vector< string > fields;
		boost::split(fields, line, boost::is_any_of("\t"));
		chrom = fields[colIdForChr];
		if (chrom != chr) continue;
		else {
			do {
				if (line.empty() || line[0]=='#') continue;
				vector< string > chrfields;
				boost::split(chrfields, line, boost::is_any_of("\t"));
				chrom = chrfields[colIdForChr];
				if (chrom != chr) break;
				start=stoi(chrfields[colIdForStart]);
				if (chrfields[colIdForCDIF]!="NA") {
					cdif=stold(chrfields[colIdForCDIF]);
					cdifs[start]=cdif;
				}
			} while ((fin.good() && !fin.eof()) && getline(fin, line));
			break;
		}
	}
	fin.close();
	return 0;
}

int readCDIF(
		string& infile
		, map< int, double >& cdifs
		, string& chr
		, int lower
		, int upper
		)
{
	int colIdForChr=0;
	int colIdForStart=1;
	int colIdForCDIF=13;

	string chrom;
	int start;
	double cdif;

	ifstream fin(infile);
	string line;
	while ((fin.good() && !fin.eof()) && getline(fin, line)) {
		if (line.empty() || line[0]=='#') continue;
		vector< string > fields;
		boost::split(fields, line, boost::is_any_of("\t"));
		chrom = fields[colIdForChr];
		if (chrom != chr) continue;
		else {
			do {
				if (line.empty() || line[0]=='#') continue;
				vector< string > chrfields;
				boost::split(chrfields, line, boost::is_any_of("\t"));
				chrom = chrfields[colIdForChr];
				if (chrom != chr) break;
				start=stoi(chrfields[colIdForStart]);
				if (start<lower) continue; // [lower, upper)
				if (start>=upper) break; // BED file should be sorted within a chromosome.
				if (chrfields[colIdForCDIF]!="NA") {
					cdif=stold(chrfields[colIdForCDIF]);
					cdifs[start]=cdif;
				}
			} while ((fin.good() && !fin.eof()) && getline(fin, line));
			break;
		}
	}
	fin.close();
	return 0;
}

int Hist::counts2hist()
{
	if (counts.empty()) return 1;
	long long N=0;
	for (long long c: counts) {
		N+=c;
	}
	if (N==0) return 1;
	for (long long c: counts) {
		pi.push_back(1.0*c/N);
	}
	double h=(upper-lower)/numbins;
	for (long long c: counts) {
		density.push_back(1.0*c/(N*h));
	}
	for (int i=0; i<numbins; ++i) {
		double left=lower+i*h;
		double right=lower+(i+1)*h;
		string label="[";
		label+=tostring(left)+",";
		label+=tostring(right);
		label+=(i==numbins-1)?"]":")";
		labels.push_back(label);
	}
	return 0;
}

int Hist::output(ostream& out)
{
	out << "label\tcounts\tdensity\tpi" << endl;
	for (int i=0; i<numbins; ++i) {
		out << labels[i] << '\t'
			<< counts[i] << '\t'
			<< density[i] << '\t'
			<< pi[i] << endl;
	}
	return 0;
}

int Hist::output(string& outfile)
{
	ofstream fout(outfile);
	output(fout);
	fout.close();
	return 0;
}

int DMR::clear()
{
	chrom.clear();
	start=-1;
	end=-1;
	meancdif.clear();
	kld.clear();
	sd.clear();
	label.clear();
	return 0;
}

bool pairgreater(pair< double, int >& a, pair< double, int >& b)
{
	return a.first > b.first;
}

int padjustfdr(vector< double >& p, vector< double >& q)
{
	vector< pair< double, int > > x;
	for (int i=0; i<p.size(); i++) {
		if (p[i]>-1) {
			x.push_back(make_pair(p[i], i));
		}
	}
	sort(x.begin(), x.end(), pairgreater);
	q.resize(p.size(), -1);
	int N=x.size();
	double minvalue=x[0].first;
	q[x[0].second]=minvalue;
	for (int i=2; i<=N; i++) {
		double value=x[i-1].first*N/(N-i+1);
		if (value<minvalue) {
			minvalue=value;
		}
		q[x[i-1].second]=(minvalue>1.0)?1.0:minvalue;
	}
	return 0;
}

