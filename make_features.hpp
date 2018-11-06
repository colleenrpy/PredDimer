#include <cstdlib>
#include <map>
#include<iostream>
#include<string>
#include<fstream>
#include <vector>
#include <set>
#include <algorithm>
#include <boost/bimap.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost::bimaps;

#define PI 3.14

struct PPIWithWeight{
  double weight;
  int node1;
  int node2;
};

struct FeatureOfWeight{
  double weight;
  double max_weight;
  double min_weight;
  double maxminw;
  double maxabsw;
  double weight_rate;
  double maxdomain;
  double mindomain;
 // double max_amino_hydroScale;
 // double maxmin_amino_hydroScale;
  double sum_weight1;
  double sum_weight2;
  double abs_max1minus2;
  double abs_min1minus2;
};

vector< FeatureOfWeight > feature;
void amino_hydrophilicity(map<char,double>& hydroScale);
void read_protnamemap(bimap<string,int>& protnamemap, vector< PPIWithWeight >& edge, vector< set< pair<int,int> > >& node, ifstream& f);
//void read_domnamemap(map< int,vector<int> >& prot_domain, map<string,int>& prot_domNum, map<string,int>& domnamemap, bimap<string,int>& protnamamap, ifstream& f);
void read_domnamemap(map< int,vector<int> >& prot_domain, map<int, vector<int> >& prot_domseq, map<int, vector<string> >&prot_location, map<int, string >& prot_sequence ,map<string,int>& prot_domNum, map<string,int>& domnamemap, bimap<string,int>& protnamamap, ifstream& f);
//void make_feature(map<string,int>& prot_domNum, bimap<string,int>& protnamemap, vector< PPIWithWeight >& edge, vector< set< pair<int,int> > >& node, vector< FeatureOfWeight>& feature);
void make_feature(map<char, double>& hydroScale, map<int,string>& prot_sequence, map<int, vector<int> >& prot_domseq, map<string,int>& prot_domNum, bimap<string,int>& protnamemap, vector< PPIWithWeight >& edge, vector< set< pair<int,int> > >& node, vector< FeatureOfWeight>& feature);

typedef bimap<string, int>::value_type value_t;

ofstream of;
//void operate(vector<PPIWithWeight>& edge, vector< set < pair<int,int> > >& node, bimap<string,int>& protnamemap, map< int,vector<int> >& prot_domain, map<string,int>& prot_domNum, map<string,int>& domnamemap) { //record proteins and indexes
void operate(vector<PPIWithWeight>& edge, vector< set < pair<int,int> > >& node, bimap<string,int>& protnamemap, map< int,vector<int> >& prot_domain, map< int,vector<int> >& prot_domseq, map<int, vector<string> >& prot_location, map<int, string>& prot_sequence, map<char,double>& hydroScale, map<string,int>& prot_domNum, map<string,int>& domnamemap) { //record proteins and indexes
#if 0
  int main(){
    vector< PPIWithWeight > edge; //record proteins and their links
    vector< set< pair<int,int> > > node;
    bimap<string,int> protnamemap;
    map<int, vector<int> >prot_domain;
    map<string,int> prot_domNum;
    map<string,int> domnamemap;
#endif
    ifstream f;
    map<string,int>::iterator it;
    map<int, vector<int> >::iterator ir;
    f.open("pro200600448_3_s.csv");
    // f.open("wi-phi-rrw-0.7");
    read_protnamemap(protnamemap, edge, node, f);
    f.close();

    f.open("uniprot_sprot_fungi.dat");
    of.open("ProtDomain");
    //read_domnamemap(prot_domain, prot_domNum, domnamemap, protnamemap, f);
    read_domnamemap(prot_domain, prot_domseq, prot_location, prot_sequence, prot_domNum, domnamemap, protnamemap, f);

    f.close();
    of.close();
    //-------------insert proteins that can't be found in file uniprot_sprot_fungi.dat by hand----
    prot_domNum.insert(pair<string,int>("YLR236C",0));
    prot_domNum.insert(pair<string,int>("YPR050C",0));
    prot_domNum.insert(pair<string,int>("YLR235C",0));
    prot_domNum.insert(pair<string,int>("YDR203W",0));
    prot_domNum.insert(pair<string,int>("YIR020C-A",0));
    prot_domNum.insert(pair<string,int>("YDR154C",0));
    prot_domNum.insert(pair<string,int>("YMR306C-A",0));
    prot_domNum.insert(pair<string,int>("YMR316C-A",0));
    prot_domNum.insert(pair<string,int>("YDR157W",0));
    prot_domNum.insert(pair<string,int>("YMR316C-B",0));
    prot_domNum.insert(pair<string,int>("YMR304C-A",0));
    prot_domNum.insert(pair<string,int>("YDR187C",0));
    prot_domNum.insert(pair<string,int>("YDR199W",0));
    prot_domNum.insert(pair<string,int>("YMR153C-A",0));
    prot_domNum.insert(pair<string,int>("YMR119W-A",0));
    prot_domNum.insert(pair<string,int>("YDR048C",0));
    prot_domNum.insert(pair<string,int>("YMR290W-A",0));
    prot_domNum.insert(pair<string,int>("YMR294W-A",0));
    prot_domNum.insert(pair<string,int>("YPR170C",0));
    prot_domNum.insert(pair<string,int>("YDR230W",0));
    prot_domNum.insert(pair<string,int>("YER187W-A",0));
    prot_domNum.insert(pair<string,int>("YMR173W-A",0));
    prot_domNum.insert(pair<string,int>("YLR458W",0));
    //-------------insert proteins that can't be found in file uniprot_sprot_fungi.dat by hand----
    //make_feature(prot_domNum, protnamemap, edge, node, feature);
    make_feature(hydroScale, prot_sequence, prot_domseq, prot_domNum, protnamemap, edge, node, feature);
}

void read_domnamemap(map< int,vector<int> >& prot_domain, map<int, vector<int> >& prot_domseq, map<int, vector<string> >& prot_location, map<int, string>& prot_sequence, map<string,int>& prot_domNum, map<string,int>& domnamemap, bimap<string,int>& protnamemap, ifstream& f){
  int domain = 0;
  int flag_of_yeast = 0;
  string line;
  string gn,former,latter;
  vector<string> prot;
  vector<int> dom;
  vector<int> domseq;
  vector<string> loca;
  string sequence;
  int flag_sequence = 0;
  while (getline(f,line) ) {
    string temp = line;
    size_t pos1,pos2,pos3;
    pos1 = temp.find(" ",0);
    if ( (temp.substr(0,pos1) == "ID")&&(temp.find("YEAST") != string::npos) ) {
      pos1 = temp.find_last_of(" ");
      pos2 = temp.find_last_of(" ", pos1-1);
      flag_of_yeast = 1;
    }
    if ( (temp.substr(0,pos1) == "GN")&&(flag_of_yeast == 1) )
    {
      pos1 = temp.find("OrderedLocusNames");
      if(pos1 != string::npos){
        pos2 = temp.find("=", pos1);
        pos3 = temp.find(";", pos2);
        gn = temp.substr(pos2+1,pos3-pos2-1);
        if (gn.find("/") != string::npos) {
          former = gn.substr(0,gn.find("/"));
          latter = gn.substr(gn.find("/")+1);
          prot.push_back(former);
          prot.push_back(latter);
        } else if(gn.find(",") != string::npos) {
          former = gn.substr(0,gn.find(","));
          latter = gn.substr(gn.find(",")+2);
          prot.push_back(former);
          prot.push_back(latter);
        } else {
          prot.push_back(gn);
        }
      }
      pos1 = temp.find("ORFNames");
      if((pos1 != string::npos)&&(flag_of_yeast == 1)){
        pos2 = temp.find("=", pos1);
        pos3 = temp.find(";", pos2);
        gn = temp.substr(pos2+1,pos3-pos2-1);
        if (gn.find("/") != string::npos) {
          former = gn.substr(0,gn.find("/"));
          latter = gn.substr(gn.find("/")+1);
          prot.push_back(former);
          prot.push_back(latter);
        } else if(gn.find(",") != string::npos) {
          former = gn.substr(0,gn.find(","));
          latter = gn.substr(gn.find(",")+2);
          prot.push_back(former);
          prot.push_back(latter);
        } else {
          prot.push_back(gn);
        }
      }
    }
    
    pos2 = temp.find(" ",pos1+5);//the first position of the second space
    if ( (temp.substr(0,pos1) == "FT")&&(temp.substr(pos1+3,pos2-pos1-3) == "DOMAIN")&&(flag_of_yeast == 1)) {
      pos3 = temp.find(" ", pos2+8);//the first position of the third space
      size_t pos6 = temp.find_last_of(" ",pos3-1);//the last position of the second space
      size_t pos4 = temp.find(" ", pos3+6); //the first position of the forth space
      size_t pos5 = temp.find_last_of(" ", pos4-1);//the last position of the thrid space
      try{
        domseq.push_back(boost::lexical_cast<int>(temp.substr(pos6+1,pos3-pos6-1)));
        domseq.push_back(boost::lexical_cast<int>(temp.substr(pos5+1,pos4-pos5-1)));
      } catch (boost::bad_lexical_cast& ){}
    }

    //domain
    pos1 = temp.find("Pfam");
    if ( (pos1 != string::npos)&&(flag_of_yeast == 1) )
    {
      pos2 = temp.find_last_of(";");
      pos3 = temp.find_last_of(".");
      size_t pos4 = temp.find(" ",9);
      size_t pos5 = temp.find(";",(int)pos4);
      string d = temp.substr((int)pos4+1,(int)pos5-(int)pos4-1);
      domnamemap.insert(pair<string,int>(d,domnamemap.size()));
      dom.push_back(domnamemap.find(d)->second);
      domain += boost::lexical_cast<int>(temp.substr((int)pos2+2,(int)pos3-(int)pos2-2));
    }

    //location
    pos1 = temp.find("CC   -!- SUBCELLULAR LOCATION:");
    if ( (pos1 != string::npos)&&(flag_of_yeast == 1) ){
      string location;
      pos2 = temp.find(":");
      pos3 = temp.find(" ",int(pos2));
      string caps(",.;");
      while((temp.find_first_of(caps, int(pos3)) != string::npos)){
        size_t pos4 = temp.find_first_of(caps,int(pos3));
        if(pos4 != string::npos) {
          location = temp.substr((int)pos3+1, (int)pos4-(int)pos3-1);
          loca.push_back(location);
        }
        pos3 = (int)pos4 + 1;
      }
    }
    //sequence
    pos1 = temp.find("SQ   SEQUENCE");
    if ( (pos1 != string::npos)&&(flag_of_yeast == 1) ){
      flag_sequence = 1;
    }
    if((flag_sequence == 1)&&(temp.find(";")==string::npos)&&(temp != "//")){
      sequence += temp;
    }

    if((temp == "//")&&(flag_of_yeast == 1)) {
      flag_sequence = 0;

      /*bubble sort*/
      if(dom.size() > 1){
        for(int j = 0; j < dom.size()-1; j++){
          for(int k = j+1; k < dom.size(); k++){
            int term;
            if(dom[j] > dom[k]){
              term = dom[j];
              dom[j] = dom[k];
              dom[k] = term;
            }
          }
        }
      }

      if(loca.size() > 1){
        for(int j = 0; j < loca.size()-1; j++){
          for(int k = j+1; k < loca.size(); k++){
            string term;
            if(loca[j] > loca[k]){
              term = loca[j];
              loca[j] = loca[k];
              loca[k] = term;
            }
          }
        }
      }

      for(int i = 0; i < prot.size(); i++){
        //       if(domain == 0) domain = 1;
        prot_domNum.insert(pair<string,int>(prot[i],domain));
        if(protnamemap.left.find(prot[i]) != protnamemap.left.end()){
          of<<prot[i]<<": ";
          int key = protnamemap.left.at(prot[i]);
          for(int j = 0; j < dom.size(); j++){
            prot_domain[key].push_back(dom[j]);
            of<<dom[j]<<" ";
          }
          of<<endl;
#if 1
          for(int j = 0; j < loca.size(); j++){
            prot_location[key].push_back(loca[j]);
            of<<loca[j]<<" ";
          }
          of<<endl;
          for(int j = 0; j < domseq.size(); j++){
            prot_domseq[key].push_back(domseq[j]);
            of<<domseq[j]<<" ";
          }
          of<<endl;
#endif
          size_t pos;
          while((pos = sequence.find_first_of(" ")) != string::npos)
          {
            sequence.erase(pos,1);
          }

          prot_sequence[key] = sequence;
          of<<sequence<<endl;
          of<<endl;
        }
      }
      dom.clear();
      loca.clear();
      sequence.clear();
      domseq.clear();
      read_domnamemap(prot_domain, prot_domseq, prot_location, prot_sequence, prot_domNum,domnamemap,protnamemap, f);
    }
  }
}
#if 1
void read_protnamemap(bimap<string,int>& protnamemap, vector< PPIWithWeight >& edge, vector< set< pair<int,int> > >& node, ifstream& f){
  size_t pos1,pos2,pos3;
  string str = "\t";	
  int edge_index = 0;
  char line[10000];
  while( f.getline(line, sizeof(line)) ){
    string temp = line;
    pos1=temp.find(str,0);
    pos2=temp.find(str,pos1+1);
    pos3=temp.find(str,pos2+1);
#if 1 //original
    if (temp.find("Score") == string::npos){
      double w;
      try{
        w = boost::lexical_cast<double>(temp.substr(0,int(pos1)));
      } catch (boost::bad_lexical_cast& ){}
      string n1str = temp.substr(int(pos1)+2,int(pos2)-3-int(pos1));
      string n2str = temp.substr(int(pos2)+2,int(pos3)-3-int(pos2));
#endif

#if 0 //modified for "wi-phi-rrw-0.7" file
      if(temp.find(str) != string::npos){
        string n1str = temp.substr(0,int(pos1));
        string n2str = temp.substr(int(pos1)+1, int(pos2)-1-int(pos1));
        double w;
        try{
          w = boost::lexical_cast<double>(temp.substr(int(pos2)+1));
        } catch (boost::bad_lexical_cast& ){}
        //   cerr<<n1str<<" "<<n2str<<" "<<w<<endl;
#endif
        if(n1str != n2str){//not consider self-loop 
          int n1;
          int n2;
          if (protnamemap.left.find(n1str) != protnamemap.left.end()) {   //if we find the name of protein1 in protnamemap
            n1 = protnamemap.left.at(n1str);
            if(protnamemap.left.find(n2str) == protnamemap.left.end()) {  //if we can't find the name of protein2 in protnamemap
              n2 = protnamemap.size();
              protnamemap.insert(value_t(n2str,n2));
              node.resize(n2+1);
              node[n2].insert(pair<int,int>(edge_index,n1));
              node[n1].insert(pair<int,int>(edge_index,n2));
            }else {
              n2 = protnamemap.left.at(n2str);
              node[n2].insert(pair<int,int>(edge_index,n1));
              node[n1].insert(pair<int,int>(edge_index,n2));
            }
          }else{  //if we can't find the name of protein1 in protnamemap
            n1 = (int)protnamemap.size();
            protnamemap.insert(value_t(n1str,n1));
            node.resize(n1+1);
            if (protnamemap.left.find(n2str) != protnamemap.left.end()) {  // if we find the name of protein2 in protnamemap
              n2 = protnamemap.left.at(n2str);
              node[n1].insert(pair<int,int>(edge_index,n2));
              node[n2].insert(pair<int,int>(edge_index,n1));
            }else{  //if we can't find the name of protein2 in protnamemap
              n2 = (int)protnamemap.size();
              protnamemap.insert(value_t(n2str,n2));
              node.resize(n2+1);
              node[n2].insert(pair<int,int>(edge_index,n1));
              node[n1].insert(pair<int,int>(edge_index,n2));
            }
          }
          PPIWithWeight ppitemp;
          ppitemp.weight = w;
          ppitemp.node1 = n1;
          ppitemp.node2 = n2;
          edge.push_back(ppitemp);
          edge_index++;
        }
      }
    }
}

void make_feature(map<char,double>& hydroScale, map<int, string>& prot_sequence, map<int,vector<int> >& prot_domseq, map<string,int>& prot_domNum, bimap<string,int>& protnamemap, vector< PPIWithWeight >& edge, vector< set< pair<int,int> > >& node, vector<FeatureOfWeight>& feature){
//void make_feature(map<string,int>& prot_domNum, bimap<string,int>& protnamemap, vector< PPIWithWeight >& edge, vector< set< pair<int,int> > >& node, vector<FeatureOfWeight>& feature){
  FeatureOfWeight fw;
  set< pair<int,int> >::iterator is;
  set< pair<int,int> >::iterator ib;
  for(int i = 0;i < edge.size();i++){
    double max_w = 0.0;
    double min_w = node[edge[i].node1].begin()->first;
    double minweight = 0.0;
    double maxweight = 0.0;
    double maxabsw = 0.0;
    double w_rate = 0.0;
    double maxdom = 0.0;
    double mindom = 0.0;
    double inner_product_maxsum1_2 = 0.0;
 //   double maxmin_amino_hydro = 0.0;
    double sum_w1 = 0.0;
    double sum_w2 = 0.0;
    double abs_max1m2 = 0.0;
    double asb_min1m2 = 0.0;

    //---------------create feature weight_rate-----------------------
    double sum_weight = 0.0;
    double max1 = 0;
    double max2 = 0;
    double min1 = 1000;
    double min2 = 1000;
    for(is = node[edge[i].node1].begin(); is != node[edge[i].node1].end(); ++is) {
      sum_weight += edge[is->first].weight;
      if(edge[is->first].weight > max1) max1 = edge[is->first].weight;
      if(edge[is->first].weight < min1) min1 = edge[is->first].weight;
    }
    sum_w1 = sum_weight;
    for(ib = node[edge[i].node2].begin(); ib != node[edge[i].node2].end(); ++ib) {
      sum_weight += edge[ib->first].weight;
      sum_w2 += edge[ib->first].weight;
      if(edge[is->first].weight > max2) max2 = edge[is->first].weight;
      if(edge[is->first].weight < min2) min2 = edge[is->first].weight;
    }
    w_rate = edge[i].weight/(sum_weight-edge[i].weight); 
    fw.weight_rate = w_rate;
    //---------------create feature weight_rate-----------------------

    //-----------create feature maxdomain and mindomain for each edge----------
    double domain1,domain2;
    string prot1 = protnamemap.right.at(edge[i].node1);
    if (prot1.find("|") != string::npos) {
      string former = prot1.substr(0,prot1.find("|"));
      string latter = prot1.substr(prot1.find("|")+1);
      if ((prot_domNum.find(former) != prot_domNum.end())&&(prot_domNum.find(latter) != prot_domNum.end())) {
        domain1 = (prot_domNum.find(former)->second + prot_domNum.find(latter)->second)/2;
      }        
    } else {
      if (prot_domNum.find(prot1) != prot_domNum.end()) {
        domain1 = prot_domNum.find(prot1)->second;
      } else {
        cout<<"Can't find the first protein "<<prot1<<" in prot_domNum"<<endl;
      }
    }
    string prot2 = protnamemap.right.at(edge[i].node2);
    if (prot2.find("|") != string::npos) {
      string former = prot2.substr(0,prot2.find("|"));
      string latter = prot2.substr(prot2.find("|")+1);
      if ((prot_domNum.find(former) != prot_domNum.end())&&(prot_domNum.find(latter) != prot_domNum.end())) {
        domain2 = (prot_domNum.find(former)->second + prot_domNum.find(latter)->second)/2;
      } 
    } else {
      if (prot_domNum.find(prot2) != prot_domNum.end()) {
        domain2 = prot_domNum.find(prot2)->second;
      } else {
        cout<<"Can't find the second protein "<<prot2<<" in prot_domNum"<<endl;
      }
    }

    if(domain1 >= domain2) {
      maxdom = domain1;
      mindom = domain2;
    } else {
      maxdom = domain2;
      mindom = domain1;
    }
    //-----------create feature maxsize_self and minsize_self for each edge----------

    for(is = node[edge[i].node1].begin(); is != node[edge[i].node1].end(); ++is) {
      //---------------create feature max{min{w1,w2}}------------------
      for(ib = node[edge[i].node2].begin(); ib != node[edge[i].node2].end(); ++ib) {
        if(((*is).second == (*ib).second)&&(edge[i].node1 != (*ib).second)&&(edge[i].node2 != (*is).second)) {
          if(edge[(*is).first].weight <= edge[(*ib).first].weight) {
            minweight = edge[(*is).first].weight;
          }else {
            minweight = edge[(*ib).first].weight;
          }
          double absw = abs(edge[(*is).first].weight - edge[(*ib).first].weight);
          if(maxabsw < absw)
            maxabsw = absw;
        }
      }
      if (minweight > maxweight) {maxweight = minweight;}
      if((*is).first != i){
        if(edge[(*is).first].weight > max_w) max_w = edge[(*is).first].weight;
        if(edge[(*is).first].weight < min_w) min_w = edge[(*is).first].weight;
      }
    }
    if(minweight == 0.0) maxweight = 0.0;
    //if(maxweight == 0.0) cout<<"There is no nodes that not only connect with node1 but also connect with node2."<<endl;
    fw.maxminw = maxweight;
    //---------------create feature max{min{w1,w2}}------------------

    for(is = node[edge[i].node2].begin(); is != node[edge[i].node2].end(); ++is) {
      if((*is).first != i){
        if(edge[(*is).first].weight > max_w) max_w = edge[(*is).first].weight;
        if(edge[(*is).first].weight < min_w) min_w = edge[(*is).first].weight;
      }
    }
#if 0
    //---------------create feature amino_hydrophilicity-------------
    int domseq_b;//beginning of domain sequence
    int domseq_e;//end of domain sequence
    string seq1,seq2;
    double scale_domseq1,scale_domseq2, scale_domseq3;
    double inner_product_sum1_2 = 0;
    if((prot_sequence.find(edge[i].node1) != prot_sequence.end())&&(prot_sequence.find(edge[i].node2) != prot_sequence.end())){
      seq1 = prot_sequence.find(edge[i].node1)->second;
      seq2 = prot_sequence.find(edge[i].node2)->second;
      int j=0;
      vector<double> sumscale_domseq1;
      vector<double> sumscale_domseq2;
      double max1 = 0;
    //  while(j < prot_domseq.find(edge[i].node1)->second.size()){
        domseq_b = (prot_domseq.find(edge[i].node1)->second)[j];
        domseq_e = (prot_domseq.find(edge[i].node1)->second)[j+1];
        scale_domseq1 = 0.0;
        int cnt1 = 1;
      //  for(int k = domseq_b; k < domseq_e+1; k++){
        for(int k = 1; k < seq1.size()+1; k++){
          double s = hydroScale.find(seq1.at(k-1))->second;
          scale_domseq1 += s*cos(PI*cnt1);
          cnt1++;
        //  if(s > 0) cnt1++;
        }
        //scale_domseq1 = (double)cnt1/(domseq_e+1-domseq_b);
        if (scale_domseq1 > max1) max1 = scale_domseq1;
       // sumscale_domseq1.push_back(scale_domseq1);
        j = j+2;
  //    }
      j = 0;
      double max2 = 0;
    //  while(j < (prot_domseq.find(edge[i].node2)->second).size()){
        domseq_b = (prot_domseq.find(edge[i].node2)->second)[j];
        domseq_e = (prot_domseq.find(edge[i].node2)->second)[j+1];
        scale_domseq2 = 0.0;
        int cnt2 = 1;
        //for(int k = domseq_b; k < domseq_e+1; k++){
        for(int k = 1; k < seq2.size()+1; k++){
          double s = hydroScale.find(seq2.at(k-1))->second;
          scale_domseq2 += s*cos(PI*cnt2);
          cnt2++;
         // scale_domseq2 += s;
        }
       // scale_domseq2 = (double)cnt2/(domseq_e+1-domseq_b);
       // sumscale_domseq2.push_back(scale_domseq2);
        if (scale_domseq2 > max2) max2 = scale_domseq2;
        j = j+2;
    //  }
      inner_product_sum1_2 = max1+max2;

      for(is = node[edge[i].node1].begin(); is != node[edge[i].node1].end(); ++is) {
        for(ib = node[edge[i].node2].begin(); ib != node[edge[i].node2].end(); ++ib) {
          if(((*is).second == (*ib).second)&&(edge[i].node1 != (*ib).second)&&(edge[i].node2 != (*is).second)) { //if have common node
            string seq3;
            double min_amino_hydro = 0.0;
            if((prot_sequence.find((*is).second) != prot_sequence.end())){
              seq3 = prot_sequence.find((*is).second)->second;
              int j = 0;
              vector<double> sumscale_domseq3;
              double max3 = 0;
             // while(j < (prot_domseq.find((*is).second)->second).size()){
                domseq_b = (prot_domseq.find((*is).second)->second)[j];
                domseq_e = (prot_domseq.find((*is).second)->second)[j+1];
                scale_domseq3 = 0.0;
                int cnt3 = 1;
                //for(int k = domseq_b; k < domseq_e+1; k++){
                for(int k = 1; k < seq3.size()+1; k++){
                  double s = hydroScale.find(seq3.at(k-1))->second;
                  scale_domseq3 += s*cos(PI*cnt3);
                  cnt3++;
              //    scale_domseq3 += s;
                }
                //scale_domseq3 = (double)cnt3/(domseq_e+1-domseq_b);
                if (scale_domseq3 > max3) max3 = scale_domseq3;
                //sumscale_domseq3.push_back(scale_domseq3);
                j = j+2;
           //   }
              double inner_product_sum1_3 = 0;
              inner_product_sum1_3 = max1+max3;
              double inner_product_sum2_3 = 0;
              inner_product_sum2_3 = max2+max3;
              if(inner_product_sum1_3 < inner_product_sum2_3){
                min_amino_hydro = inner_product_sum1_3;//sum_high_hydro2/cnt2;
              }else{
                min_amino_hydro = inner_product_sum2_3;//sum_high_hydro1/cnt1;
              }
            } else{
              cerr<<"seq3 is not exist. seq3 = "<<seq3<<endl;
            }
            if(min_amino_hydro > maxmin_amino_hydro) maxmin_amino_hydro = min_amino_hydro;
          }
        }
      }
    }
    // if(ave_seq1 > ave_seq2) amino_hydro = ave_seq1;
    // else amino_hydro = ave_seq2;
#endif
    fw.weight = edge[i].weight;
    fw.max_weight = max_w;
    fw.min_weight = min_w;
    fw.maxabsw = maxabsw;
    fw.maxdomain = maxdom;
    fw.mindomain = mindom;
  //  fw.max_amino_hydroScale = inner_product_sum1_2;
  //  fw.maxmin_amino_hydroScale = maxmin_amino_hydro;
    fw.sum_weight1 = sum_w1;
    fw.sum_weight2 = sum_w2;
    fw.abs_max1minus2 = abs(max1-max2);
    fw.abs_min1minus2 = abs(min1-min2);
    feature.push_back(fw);
  }
}
#endif
