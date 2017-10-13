#include"pFind_PairResearch.h"
#include<algorithm>
#include<unordered_map>
#include<fstream>
#include<iostream>
static const int Utf8_ONE = 49;//UTF-8�У��ַ���1���ı�����49
ostream & operator<<(ostream & os, const spectra & s)
{
	os << s.file_name << "	" << s.scan_no << "	" << s.exp_mh << "	" << s.charge << "	" << s.q_value << "	" << s.seq << "	"
		<< s.calc_mh << "	" << s.mass_shift << "	" << s.raw_score << "	" << s.final_score << "	" << s.modi << "	" << s.spec
		<< "	" << s.prot << "	" << s.posi << "	" << s.label << "	" << s.targe << "	" << s.mc_sites << "	" << s.afm_shift << "	"
		<< s.others;
	//����������ͻ����ʱ�����ð������
	os <<"	"<<s.marker<< "	"<<s.is_mut;//��������ͻ���ļ��У�ͻ����Ϣ������ENSP��ĩβ���ֶΣ�_SAP),��ȡensp��ʱ���е����ֶΣ�Ϊ������ļ�������ͻ����ͬ���������������һ���������
	return os;
}

istream & operator>>(istream & in, vector<spectra> & list_result)
{
	char ch;
	string waste;
	getline(in, waste);
	spectra temp;//��Ȼ��ÿ�γ�ʼ���������Ż��»���ÿ�γ�ʼ����ռ��ͬһ�ڴ����򣬴Ӷ�ʹis_mut,modi�ȱ���
	//ֱ����ֵǰ���ᱣ����һ�ε�ֵ����Ϊ�ڴ�ǡ�ö����ˣ�����˶Ժ��������о��������õ�is_mutֵ����Ӧ����������if else������¸�ֵ
	while (in >> temp.file_name) {
		in >> temp.scan_no;
		in >> temp.exp_mh;
		in >> temp.charge;
		in >> temp.q_value;
		in >> temp.seq;
		in >> temp.calc_mh;
		in >> temp.mass_shift;
		in >> temp.raw_score;
		in >> temp.final_score;
		in.get(ch).get(ch);
		if (ch == '\t')
			temp.modi = "";
		else
		{
			in >> temp.modi;
			//����������ͻ����ʱ�����ð������
			temp.modi = ch + temp.modi;//��openresearch�ļ��У���������һ��˫���Ű�����Ϊ��֤��ȡͻ����Ϣʱ���㣬�˴�����ǰ˫���ţ�������ͻ�����������û��˫����
			//������openresearch���ʱ�����ð������
			///if (temp.modi.find(">") != string::npos) {
			///	temp.is_mut = true;
			///}
			///else
			///	temp.is_mut = false;
		}
		in >> temp.spec;
		in >> temp.prot;
		//����������ͻ����ʱ�����ð������
		if (temp.prot.find("_SAP") != string::npos)
		{
			temp.is_mut = true;
			temp.prot.erase(temp.prot.find("_SAP"), 4);
		}
		else
			temp.is_mut = false;
		temp.prot.erase(temp.prot.length() - 1,1);//�������һλ�ķ�б��
		in >> temp.posi;
		in >> temp.label;
		in >> temp.targe;
		in >> temp.mc_sites;
		in >> temp.afm_shift;
		in >> temp.others;
		//������openresearch���ʱ�����ð������
		///in >> temp.mut_count;
		if(temp.prot.length()<25&&temp.targe=="target"&&temp.q_value<0.01)
			list_result.push_back(temp);
	}
	return in;
}

int mark(vector<spectra>& list)//�����ݷֳɼ������飬ÿ������ͬ����С�����У��Ҳ��ظ������˺���δ���ͬʱ����������ͬ����С�����У��Ӷ����Թ����ڶ����ͬ�飩�����⡣��Ը���ݼ���û���������Ķ�����
{
	sort(list.begin(), list.end(), sortbyleg);
	int mark = 1;//ע�⣬mark��1��ʼ��
	for (int i = 0; i < list.size(); i++) {
		if (list[i].marker != 0)
			continue;
		else
			list[i].marker = mark++;//��ֵ�����
		for (int j = i + 1; j < list.size(); j++) {
			if (list[j].marker != 0)
				continue;
			if (list[j].seq.find(list[i].seq) != string::npos)//����ҵ��˵Ļ�
				list[j].marker = list[i].marker;
		}
	}
	sort(list.begin(), list.end(), sortbymarker);
	return mark;//����mark�����ڽ���bool�����ж��Ƿ�������
}



bool sortbyleg(spectra & a, spectra & b)
{
	return a.seq.length()<b.seq.length();
}

bool sortbymarker(spectra & a, spectra & b)
{
	return a.marker<b.marker;
}

mut_pep_inform pepmutation(const spectra & p, ifstream& intri, unordered_map<string, char>& table)
{
	//�ı�����
	string modi = p.modi;
	mut_pep_inform a;
		a.mutpep = p.seq;
		a.size = p.mut_count;
		a.pos_mut = new int[a.size]{ 0 }; 
	string section;
	int order = 0;
	while (modi.find(">") != string::npos && modi.find(";") != string::npos) {
		int end_sign = modi.find(";");
		section = modi.substr(0, end_sign + 1);
		modi.erase(0, end_sign + 1);
		for (int i = 0; section[i] != ','; i++)//�õ�ͻ�����꣬����pos_mut��
			a.pos_mut[order] = a.pos_mut[order] * 10 + section[i] - Utf8_ONE;//���ַ������ȥ��1���ı��룬�õ��Ĳ�ֵ���Ǵ�0��ʼ�����ͻ��λ�����确8��(��1������ַ���ͻ��λ)-��1��=7(��0������ַ���ͻ��λ)
		section.erase(0, section.find(',') + 1);
		string ori_res, mut_res;
		ori_res=section.substr( 0, 3);
		int mut_sign = section.find(">");
		mut_res=section.substr(mut_sign + 1, 3); 
		char ori, mut;
		if (table.find(ori_res) != table.end() && table.find(mut_res) != table.end()) {
			ori = table[ori_res];
			mut = table[mut_res];
		}
		else
			std::cout << "���Ҳ�����д�Ĳл�������Ϊ��" <<mut_res << std::endl;
		if (a.mutpep[a.pos_mut[order]] == ori)
			a.mutpep[a.pos_mut[order]] = mut;
		order++;
		ori_res.~basic_string();
		mut_res.~basic_string();
	}
	section.~basic_string();
	return a;
}

struct tri_hash{
	size_t operator()(const string & str)const {
		int _h = 0;
		for (int i = 0; i < str.length(); i++) {
			_h = _h * 2 + str[i] - 65;
		}
		return _h;
	}
};

struct tri_compare{
	bool operator()(const string& a1, const string& a2)const {
		return a1 == a2;
	}
};