#pragma once
#include<iostream>
#include<vector>
#include<string>
#include<fstream>
const int MAX = 10;
using namespace std;
struct pro {
	std::string ensp="";
	char origaa;
	char mutataa;
	int pos;
	std::string hseq="";
};

istream& operator>>(istream& is, vector<pro>& list);
bool sortbyensp(pro& a, pro& b);
void fillhseq(istream& is,vector<pro>& list);