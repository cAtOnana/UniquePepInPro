#include "reflect.h"
#include<algorithm>
#include<unordered_map>
struct str_hash {//自定义hash函数
	size_t operator() (const string& str) const {
		unsigned long _h = 0;
		for (int i = 0; i < str.size(); i++)
			_h = _h * 2 + str[i] - 49;
		return size_t(_h);
	}
};
struct str_compare {
	bool operator()(const string& a1, const string& a2) const {
		return a1 == a2;
	}
};
typedef unordered_map<string, string> NmEnspMap,EnspHseqMap;
//////////////////////////////////////////////////////////////
istream & operator>>(istream & is, vector<pro>& list)
{
	pro temp;
	string waste;
	getline(is, waste);
	getline(is, waste);
	while (is >> waste) {
		for (int i = 0; i < 12; i++)
			is >> waste;
		is >> temp.ensp;
		is >> temp.origaa;
		is >> temp.mutataa;
		is >> temp.pos;
		getline(is,waste);
		list.push_back(temp);
	}
	return is;
}
bool sortbyensp(pro & a, pro & b)
{
	return a.ensp < b.ensp;
}

void fillhseq(istream & is, vector<pro>& list)
{
	EnspHseqMap ehmap;
	string ensp, hseq;
	while (is) {//填表
		while (is.get() != '>'&&is)//找'>'
			continue;
		is >> ensp >> hseq;
		ehmap[ensp] = hseq;
	}
	for (int i = 0; i < list.size(); i++) {
		if (ehmap.find(list[i].ensp) != ehmap.end())
			list[i].hseq = ehmap[list[i].ensp];
	}
}
