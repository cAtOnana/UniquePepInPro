#include"reflect.h"
#include"pFind_PairResearch.h"
#include<string>
#include<fstream>
#include<vector>
#include<iostream>
#include<unordered_map>

using namespace std;
string proname{ "非同义突变.txt" };
string spename{ "MHCCLM3_1.spectra" };
string outname{ "MHCCLM3_1.uniquestatic" };
int main() {
	ifstream inpro(proname);
	vector<pro> prolist;
	inpro >> prolist;

	ifstream inspe(spename);
	vector<spectra> spelist;
	inspe>>spelist;

	unordered_map<string, bool> hash_pro;
	for (auto n : prolist)
	{
		hash_pro[n.ensp] = true;
	}

	ofstream out(outname);
	for (auto n : spelist)
	{
		if (hash_pro.find(n.prot) != hash_pro.end())
			out << n << endl;
	}

	return 0;
}