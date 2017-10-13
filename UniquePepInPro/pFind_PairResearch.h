#pragma once
#include<string>
#include<vector>
#include<unordered_map>
using namespace std;
struct spectra
{
	string file_name;
	int scan_no;
	double exp_mh;
	int charge;
	double q_value;
	string seq;
	double calc_mh;
	double mass_shift;
	double raw_score;
	string final_score;
	string modi;
	int spec;
	string prot;
	string posi;
	string label;
	string targe;
	int mc_sites;
	double afm_shift;
	int others;
	int marker = 0;
	bool is_mut = false;
	bool outputable = true;
	int mut_count=0;//记录突变位点个数
};
struct mut_pep_inform {
	string mutpep;
	int* pos_mut = nullptr;
	int size=0;
};
ostream& operator<<(ostream& os, const spectra& s);
istream& operator>>(istream& os, vector<spectra>& list_result);
int mark(vector<spectra>& list);
bool sortbyleg(spectra& a, spectra& b);
bool sortbymarker(spectra& a, spectra& b);
mut_pep_inform pepmutation(const spectra& p,ifstream& intri, unordered_map<string, char>& table);
