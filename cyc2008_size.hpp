#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <boost/bimap/bimap.hpp>
#include "make_features.hpp"


using namespace std;
using namespace boost::bimaps;

typedef bimap<string,int>::value_type value_t;

void cyc2008_size(vector<vector<int> >& vector_2_size2, vector<vector<int> >& vector_2_size3, bimap<string,int>& protnamemap, map< int,vector<int> >& prot_domseq, map<int, vector<string> >& prot_location, map<int, string>& prot_sequence, map<char,double>& hydroScale, map< int, vector<int> >& prot_domain, vector<PPIWithWeight>& edge, vector<set<pair<int,int> > >& node) {
	string file1 = "cyc2008.dat";
	ifstream f;
	char line[10000];
	vector<vector<string> > cyc2008;//i record the ith name of complex and then the names of its proteins
	vector<int> Num;//i+1 record the size of ith complex
  vector<string> prot;//i record the names of proteins of ith complex for a while
  map< string, int > prot_domNum;
  map< string, int > domnamemap;
#if 0
  hydroScale.insert(make_pair('A',0.78));
  hydroScale.insert(make_pair('L',0.56));
  hydroScale.insert(make_pair('R',1.58));
  hydroScale.insert(make_pair('K',1.10));
  hydroScale.insert(make_pair('N',1.20));
  hydroScale.insert(make_pair('M',0.66));
  hydroScale.insert(make_pair('D',1.35));
  hydroScale.insert(make_pair('F',0.47));
  hydroScale.insert(make_pair('C',0.55));
  hydroScale.insert(make_pair('P',0.69));
  hydroScale.insert(make_pair('Q',1.19));
  hydroScale.insert(make_pair('S',1.00));
  hydroScale.insert(make_pair('E',1.45));
  hydroScale.insert(make_pair('T',1.05));
  hydroScale.insert(make_pair('G',0.68));
  hydroScale.insert(make_pair('W',0.70));
  hydroScale.insert(make_pair('H',0.99));
  hydroScale.insert(make_pair('Y',1.00));
  hydroScale.insert(make_pair('I',0.47));
  hydroScale.insert(make_pair('V',0.51));

  hydroScale.insert(make_pair('A',-0.5));
  hydroScale.insert(make_pair('L',-1.8));
  hydroScale.insert(make_pair('R',3.0));
  hydroScale.insert(make_pair('K',3.0));
  hydroScale.insert(make_pair('N',0.2));
  hydroScale.insert(make_pair('M',-1.3));
  hydroScale.insert(make_pair('D',3.0));
  hydroScale.insert(make_pair('F',-2.5));
  hydroScale.insert(make_pair('C',-1.0));
  hydroScale.insert(make_pair('P',0.0));
  hydroScale.insert(make_pair('Q',0.2));
  hydroScale.insert(make_pair('S',0.3));
  hydroScale.insert(make_pair('E',3.0));
  hydroScale.insert(make_pair('T',-0.4));
  hydroScale.insert(make_pair('G',0.0));
  hydroScale.insert(make_pair('W',-3.4));
  hydroScale.insert(make_pair('H',-0.5));
  hydroScale.insert(make_pair('Y',-2.3));
  hydroScale.insert(make_pair('I',-1.8));
  hydroScale.insert(make_pair('V',-1.5));
#endif
  operate(edge,node,protnamemap,prot_domain,prot_domseq, prot_location, prot_sequence, hydroScale, prot_domNum,domnamemap);

	f.open(file1.c_str());
	string str="	";
	size_t pos1,pos2,pos3;
	int cnt;
	while( f.getline(line, sizeof(line)) ) {
		string temp = line;
		pos1=temp.find(str,1);
		pos2=temp.find(str,pos1+1);
		pos3=temp.find(str,pos2+1);
		string CompName;
		if (pos3 != string::npos){
      if(prot.size()!=0){
        cyc2008.push_back(prot);
        prot.clear();
      }
			for(int i=pos2+1;i<pos3;i++){
				CompName+=temp.at(i);
			}
			prot.push_back(CompName);
			Num.push_back(cnt);
			cnt=1;
		}
		else cnt++;
    string ProtName;
		if (pos1!=string::npos){
			for(int i=0;i<pos1;i++){
				ProtName+=temp.at(i);
			}
			prot.push_back(ProtName);
		}
	}
  cyc2008.push_back(prot);
  prot.clear();
	f.close();
	Num.push_back(cnt);
  for(int i = 0; i < cyc2008.size(); i++) {
    for(int j = 1; j < cyc2008[i].size()-1; j++) {
      for(int k = j+1; k < cyc2008[i].size(); k++) {
        if (cyc2008[i][j] > cyc2008[i][k]) {
          string t;
          t = cyc2008[i][j];
          cyc2008[i][j] = cyc2008[i][k];
          cyc2008[i][k] = t;
        }
      }
    }
  }
//  typedef pair <string, int> Pair;

  //---------deal with complexes of size 2----------
  vector<int> vector_1_size2;
  bimap<string,int>::iterator p1,p2,p3;
  int num = protnamemap.size();
	ofstream of;
  of.open("CompOfSize2.txt");
  for(int i = 0; i<Num.size(); i++){
      if(Num[i] == 2) {
        of<<cyc2008[i-1][1]<<" "<<cyc2008[i-1][2]<<endl;
        if (protnamemap.left.find(cyc2008[i-1][1]) == protnamemap.left.end()) {
          protnamemap.insert(value_t(cyc2008[i-1][1],num));
          vector_1_size2.push_back(num);
          num++;
        }else {
          vector_1_size2.push_back(protnamemap.left.at(cyc2008[i-1][1]));
        }
        if (protnamemap.left.find(cyc2008[i-1][2]) == protnamemap.left.end()) {
          protnamemap.insert(value_t(cyc2008[i-1][2],num));
          vector_1_size2.push_back(num);
          num++;
        }else{
          vector_1_size2.push_back(protnamemap.left.at(cyc2008[i-1][2]));
        }
        vector_2_size2.push_back(vector_1_size2);
        vector_1_size2.clear();
      }
  }
  of.close();  
  //---------deal with complexes of size 2----------

  //---------deal with complexes of size lager than 2--------
  vector<int> vector_1_size3;
  num = protnamemap.size();
  of.open("CompOfSizeLagerThan3.txt");
  for (int i = 3; i < 82; i++) {
    for(int j = 0; j<Num.size(); j++) {
      if(Num[j] == i) {
        for(int k = 1; k <= i; k++) {
          of<<cyc2008[j-1][k]<<" ";
          if (protnamemap.left.find(cyc2008[j-1][k]) == protnamemap.left.end()) {
            protnamemap.insert(value_t(cyc2008[j-1][k],num));
            vector_1_size3.push_back(num);
            num++;
          }else {
            vector_1_size3.push_back(protnamemap.left.at(cyc2008[j-1][k]));
          }
        }
        of<<endl;
        vector_2_size3.push_back(vector_1_size3);
        vector_1_size3.clear();
      }
    }
  }
  of.close();
  //---------deal with complexes of size lager than 2--------
}




