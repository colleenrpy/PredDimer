#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <boost/bimap/bimap.hpp>

using namespace std;
using namespace boost::bimaps;

typedef bimap<string, int>::value_type value_t;

void myfeature( map< pair< string,string >, int >& compnamemap){
 // map<pair<string,string>, int> compnamemap;
  ifstream f;
  ofstream of;
  char line[10000];
  f.open("pro200600448_3_s.csv");
//  of.open("OutputMyFeature.txt");
  int i = 0;
  while(f.getline(line, sizeof(line))){
    string temp = line;
    size_t pos1,pos2,pos3;
    string str = "\t";	
    pos1=temp.find(str,0);
    pos2=temp.find(str,pos1+1);
    pos3=temp.find(str,pos2+1);
    string n1str,n2str;
    if (temp.find("Score") == string::npos){
      n1str = temp.substr(int(pos1)+2,int(pos2)-3-int(pos1));
      n2str = temp.substr(int(pos2)+2,int(pos3)-3-int(pos2));
      if(n1str != n2str){
        if(n1str<n2str) compnamemap.insert(make_pair(make_pair(n1str,n2str),i));
        else compnamemap.insert(make_pair(make_pair(n2str,n1str),i));
        i++;
      }
    }
  }
  map<pair<string,string>, int>:: iterator iter;
/*  for(iter = compnamemap.begin(); iter != compnamemap.end(); iter++) {
    of<<(*iter).first.first<<"	"<<(*iter).first.second<<"	"<<feature[(*iter).second].weight<<"	"<<feature[(*iter).second].max_weight<<"	"<<feature[(*iter).second].min_weight<<"	"<<feature[(*iter).second].sum_weight<<"	"<<feature[(*iter).second].ave_weight<<endl;
  }*/
//  of.close();
  f.close();
}
