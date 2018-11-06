#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <boost/bimap/bimap.hpp>
#include <boost/algorithm/string.hpp>
#include "MyFeature.hpp"
#include "cyc2008_size.hpp"
//#include "Phylogenetic-profile.hpp"

using namespace std;
using namespace boost::bimaps;
using namespace boost::algorithm;

void normalize_and_add_id(const string& filename1, const string& filename2, vector<vector<double> >& pos_neg, vector<double>& diag){
  ifstream fin;
  fin.open(filename1.c_str());
  ofstream fout;
  fout.open(filename2.c_str());
  if(!fin) fout<<"The file does not exist."<<endl;
  int row_id = 1;
  double normalized_value;
  int sig = pos_neg[0].size()-3;
  while((!fin.eof())&&(row_id <= pos_neg.size())){
    string line;
    vector<string> row;
    getline(fin,line);
    if(!line.empty()){
      split(row, line, is_space(), token_compress_on);
    }
#if 1
    if(pos_neg[row_id-1][sig] == 0){
      fout<<"+1 0:"<<row_id;
      for(int m = 1; m < row.size(); m++){
        //       if((diag[m-1]!=0)&&(diag[row_id-1]!=0)){
        //         normalized_value = boost::lexical_cast<double>(row[m-1])/(sqrt(diag[m-1]*diag[row_id-1]));
        //       }  //normalize
        //       else {
        normalized_value = boost::lexical_cast<double>(row[m-1]);
        //       }
        fout<<" "<<m<<":"<<normalized_value;
        //      fout<<" "<<m<<":"<<row[m-1];
      }
    }
#endif
    if(pos_neg[row_id-1][sig] == 1){
      fout<<"-1 0:"<<row_id;
      for(int m = 1; m < row.size(); m++){
        //       if((diag[m-1]!=0)&&(diag[row_id-1]!=0)){
        //         normalized_value = boost::lexical_cast<double>(row[m-1])/(sqrt(diag[m-1]*diag[row_id-1]));
        //       }  //normalize
        //       else {
        normalized_value = boost::lexical_cast<double>(row[m-1]);
        //       }
        fout<<" "<<m<<":"<<normalized_value;
        //    fout<<" "<<m<<":"<<row[m-1];
      }
    }
    fout<<endl;
    row_id++;
  } 
  fout.close();
  fin.close();
}

int main(int argc, char**argv){ 
  int feature_num = argc-4; //the number of feature
#if 1
  if (argc < 2) {
    cerr << argv[0] <<endl;
    cerr<<"features: "<<" <a>"<<" <b>"<<" <c>"<<" <d>"<<endl;
    cerr<<"coefficient: "<<argv[feature_num+1]<<endl;
    cerr<<"inputfile: "<<argv[feature_num+2]<<endl;
    cerr<<"outputfile: "<<argv[feature_num+3]<<endl;
    return 0;
  }
#endif
  //Step1:make positives and negatives
  vector<vector<int> > vector_2_size2;
  vector<vector<int> > vector_2_size3;
  bimap<string,int> protnamemap;
  vector< vector<double> > pos_neg;
  map< int, vector<int> > prot_domain;
  map< int, vector<int> > prot_domseq;
  map< int, vector<string> > prot_location;
  map< int, string> prot_sequence;
  map<pair<string,string>,int> compnamemap;
  map<char,double> hydroScale;
  map<pair<string,string>,int>::iterator cnm;
  map<pair<string,string>,int>::iterator cnm2;
  vector<PPIWithWeight> edge;
  vector< set< pair<int,int> > > node;
  cyc2008_size(vector_2_size2, vector_2_size3, protnamemap, prot_domseq, prot_location, prot_sequence, hydroScale, prot_domain, edge, node);

  myfeature(compnamemap);
  vector< vector<double> > train_train;
  vector< double > train;
  ifstream f;
  ofstream of;
  char line[10000];
  f.open("FeaturePPIWeightDiffFromMax.txt");
  size_t pos1,pos2,pos3;
  string str = "	";	
  while( f.getline(line, sizeof(line)) ){
    string temp1 = line;
    pos1 = temp1.find(str);
    pos2 = temp1.find(str,pos1+1);
    string feature = temp1.substr(pos2+1,temp1.size()-1-pos2);
    try{
      train.push_back(boost::lexical_cast<double>(feature));
    } catch(boost::bad_lexical_cast& ){}
    train_train.push_back(train);
    train.clear();
  }
  f.close();
  int i = 0; 
#if 1
  for(cnm = compnamemap.begin(); cnm != compnamemap.end(); cnm++) {
    double w_minus_min = feature[(*cnm).second].weight-feature[(*cnm).second].min_weight;
    train_train[i].push_back(feature[(*cnm).second].weight); // feature 1
    train_train[i].push_back(feature[(*cnm).second].max_weight);
    train_train[i].push_back(feature[(*cnm).second].min_weight);
    train_train[i].push_back(feature[(*cnm).second].maxminw);
    train_train[i].push_back(feature[(*cnm).second].maxabsw);
    train_train[i].push_back(feature[(*cnm).second].weight_rate);
    train_train[i].push_back(w_minus_min);
    train_train[i].push_back(feature[(*cnm).second].maxdomain);
    train_train[i].push_back(feature[(*cnm).second].mindomain);
  //  train_train[i].push_back(feature[(*cnm).second].max_amino_hydroScale);
  //  train_train[i].push_back(feature[(*cnm).second].maxmin_amino_hydroScale);
  //  train_train[i].push_back(feature[(*cnm).second].sum_weight1);
  //  train_train[i].push_back(feature[(*cnm).second].sum_weight2);
  //  train_train[i].push_back(feature[(*cnm).second].abs_max1minus2);
  //  train_train[i].push_back(feature[(*cnm).second].abs_min1minus2);
    i++;
  }
#if 1
  vector<double> tempo;

  int flag = 0;
  map<pair<string,string>, int> pos_example;
  map<pair<string,string>, int> neg_example;
  int index1 = 0;
  int index2 = 0;
  for(cnm = compnamemap.begin(); cnm != compnamemap.end(); cnm++) {
    int n1 = protnamemap.left.at((*cnm).first.first);
    int n2 = protnamemap.left.at((*cnm).first.second);
    int judge_pos=0;
    for (int j = 0; j< vector_2_size3.size(); j++){ 
      for(int k = 0; k < vector_2_size3[j].size()-1; k++) {
        for(int l = k+1; l < vector_2_size3[j].size(); l++) {
          if((n1 == vector_2_size3[j][k])&&(n2 == vector_2_size3[j][l])){ //belong to size3
            int judge_neg = 0;
            for( int i = 0; i < vector_2_size2.size(); i++){
              if((n1 == vector_2_size2[i][0])&&(n2 == vector_2_size2[i][1])) { //belong to size2
                goto H;
              } else {
                judge_neg++;
              }
            }
            if (judge_neg == 172){ // not belong to size2
              train_train[flag].push_back(1);
              train_train[flag].push_back(n1);
              train_train[flag].push_back(n2);
              if (neg_example.find(make_pair((*cnm).first.first,(*cnm).first.second)) == neg_example.end()) {
                for(int s = 0; s < train_train[flag].size(); s++){
                  tempo.push_back(train_train[flag][s]);
                }
                pos_neg.push_back(tempo);
                tempo.clear();
              }
              neg_example.insert(make_pair(make_pair((*cnm).first.first,(*cnm).first.second),index1));
              index1++;
            }
          }
          else judge_pos++;
        }
      }
    }     
    if(judge_pos == 11751) {
      for( int i = 0; i < vector_2_size2.size(); i++){
        if((n1 == vector_2_size2[i][0])&&(n2 == vector_2_size2[i][1])) {
          train_train[flag].push_back(0);
          train_train[flag].push_back(n1);
          train_train[flag].push_back(n2);
          if (pos_example.find(make_pair((*cnm).first.first,(*cnm).first.second)) == pos_example.end()) {
            for(int j = 0; j < train_train[flag].size(); j++){
              tempo.push_back(train_train[flag][j]);
            }
            pos_neg.push_back(tempo);
            tempo.clear();
          }
          pos_example.insert(make_pair(make_pair((*cnm).first.first, (*cnm).first.second),index2));
          index2++;
        }
      }
    }
H:{}
  flag++;
  }
#endif

  of.open("pos");
  for(cnm = pos_example.begin(); cnm != pos_example.end(); cnm++) {
    of<<(*cnm).first.first<<" "<<(*cnm).first.second<<endl;
  }
  of.close();
  of.open("neg");
  for(cnm = neg_example.begin(); cnm != neg_example.end(); cnm++) {
    of<<(*cnm).first.first<<" "<<(*cnm).first.second<<endl;
  }
  of.close();
  cerr<<argv[feature_num]<<endl;

  //Step2:kernels
#if 1
  double coef = boost::lexical_cast<double>(argv[feature_num+1]);
  //vector<int> feature(feature_num);
  vector<int> feature;
  cerr<<coef<<endl;
  for (int i=0; i != feature_num; ++i) {
// feature[i] = boost::lexical_cast<int>(argv[i+1]);
    feature.push_back(boost::lexical_cast<int>(argv[i+1]));
  }
  vector<double> element_row;
  vector<double> element_diagonal;
  vector< vector<double> > element_column;
  of.open(argv[feature_num+2]);
  string inputfile = "oc";
  string outputfile = "phylogenetic";
  bimap<string, double> organism_list;
  vector< vector<string> > organisms_gene;
  map<int ,vector<double> > protein_organism;
 // make_organism_map(inputfile, organism_list);
 // select_genes_organisms_with_sce(inputfile,outputfile,organism_list,organisms_gene);
 // organism_present_for_protein(protnamemap, protein_organism, organisms_gene, organism_list);

  cerr<<feature_num<<endl;
#endif

  vector<double> diag;
  int pos_prot1 = pos_neg[0].size()-2;
  int pos_prot2 = pos_neg[0].size()-1;
#if 0
  vector<double> mean_pos_f;
  vector<double> stand_pos_f;
  vector<double> mean_neg_f;
  vector<double> stand_neg_f;
  for(int i = 0; i < feature_num; i++){
    mean_pos_f.push_back(0);
    stand_pos_f.push_back(0);
    mean_neg_f.push_back(0);
    stand_neg_f.push_back(0);
  }
  //calculate mean
  for(int i = 0; i < pos_neg.size(); i++){
    for(int j = 0; j != feature_num; ++j){
      if(pos_neg[i][pos_neg[i].size()-3] == 1){
        mean_pos_f[j]+=pos_neg[i][feature[j]];
      } else {
        mean_neg_f[j]+=pos_neg[i][feature[j]];
      }
    }
  }
  for(int i = 0; i != feature_num; ++i){
    mean_pos_f[i] = mean_pos_f[i]/pos_example.size();
    mean_neg_f[i] = mean_neg_f[i]/neg_example.size();
  }
  //calculate variance
  double temp = 0;
  for(int i = 0; i < pos_neg.size(); i++){
    for(int j = 0; j != feature_num; j++){
      if(pos_neg[i][pos_neg[i].size()-3] == 1){
        stand_pos_f[j] += pow(pos_neg[i][feature[j]]-mean_pos_f[j],2);
      } else{
        stand_neg_f[j] += pow(pos_neg[i][feature[j]]-mean_neg_f[j],2);
      }
    }
  }
//  for(int i = 0; i != feature_num; ++i){
//    stand_pos_fpi] = sqrt(stand_pos_f[i]/pos_example.size());
//    stand_neg_f[i] = sqrt(stand_neg_f[i]/neg_example.size());
//  }
#endif
  ofstream ff,fff;
  fff.open("features");
  ff.open("label");
  int cnt = 0;
  for(int i = 0; i < pos_neg.size(); i++){
    if(pos_neg[i][pos_neg[i].size()-3] == 0){
      ff<<"0";
    }else{
      ff<<"1";
    }
    for(int j = 0; j < feature_num; j++){
        //ff<<" "<<(pos_neg[i][feature[j]]-mean_pos_f[j])/stand_pos_f[j];
        fff<<pos_neg[i][feature[j]];
        if(j != feature_num-1) fff<<" ";
    }
    fff<<endl;
    ff<<endl;
  }
  ff.close();
  fff.close();

  ff.open("sequence");
  ofstream fin;
  fin.open("spectrum-kernel");
  for(int i = 0; i < pos_neg.size(); i++) {
    for(int j = 0; j < pos_neg.size(); j++) {
      double mlpk_kernel = 0.0;
      double tppk_kernel = 0.0;
      double domain_kernel = 0.0;
      double element = 0.0;
#if 1
      for (int k = 0; k != feature_num; ++k) {
        element += pos_neg[i][feature[k]]*pos_neg[j][feature[k]];
      }
#endif
#if 0
      vector<double> mean_i;
      vector<double> mean_j;
      vector<double> stand_i;
      vector<double> stand_j;
      int m = pos_neg[i].size()-3;
      if(pos_neg[i][m] == 1){
        for(int k = 0; k < feature_num; k++){
          mean_i.push_back(mean_pos_f[k]);
          stand_i.push_back(stand_pos_f[k]);
        }
      } else {
        for(int k = 0; k < feature_num; k++){
          mean_i.push_back(mean_neg_f[k]);
          stand_i.push_back(stand_neg_f[k]);
        }
      }
      if(pos_neg[j][m] == 1){
        for(int k = 0; k < feature_num; ++k){
          mean_j.push_back(mean_pos_f[k]);
          stand_j.push_back(stand_pos_f[k]);
        }
      } else {
        for(int k = 0; k < feature_num; ++k){
          mean_j.push_back(mean_neg_f[k]);
          stand_j.push_back(stand_neg_f[k]);
        }
      }

      for (int k = 0; k != feature_num; ++k) {
        element += ((pos_neg[i][feature[k]]-mean_i[k])/stand_i[k])*((pos_neg[j][feature[k]]-mean_j[k])/stand_j[k]);
      }
#endif
      //comments from reviewer
      //element *= (1-coef);

#if 1
      vector<int> compare1(10000);
      vector<int> compare2(10000);
      vector<int> compare3(10000);
      vector<int> compare4(10000);
      vector<int>::iterator c1,c2,c3,c4;
      double prot1_prot3,prot1_prot4,prot2_prot3,prot2_prot4;
#if 1
      if((prot_domain.find(pos_neg[i][pos_prot1]) != prot_domain.end())&&(prot_domain.find(pos_neg[j][pos_prot1]) != prot_domain.end()))
      {
        c1 = set_intersection((prot_domain.find(pos_neg[i][pos_prot1])->second).begin(), (prot_domain.find(pos_neg[i][pos_prot1])->second).end(),(prot_domain.find(pos_neg[j][pos_prot1])->second).begin(), (prot_domain.find(pos_neg[j][pos_prot1])->second).end(),compare1.begin());
        prot1_prot3 = double(c1-compare1.begin())/((prot_domain.find(pos_neg[i][pos_prot1])->second).size()+(prot_domain.find(pos_neg[j][pos_prot1])->second).size()-double(c1-compare1.begin())); //tanimoto
        //prot1_prot3 = double(c1-compare1.begin())/sqrt((prot_domain.find(pos_neg[i][pos_prot1])->second).size()*(prot_domain.find(pos_neg[j][pos_prot1])->second).size()); //scaled
        //   prot1_prot3 = numerator1/(2-numerator1);
      }
      else {
        c1 = compare1.begin();
        prot1_prot3 = 0;
      }
      if((prot_domain.find(pos_neg[i][pos_prot1]) != prot_domain.end())&&(prot_domain.find(pos_neg[j][pos_prot2]) != prot_domain.end()))
      {
        c2 = set_intersection((prot_domain.find(pos_neg[i][pos_prot1])->second).begin(), (prot_domain.find(pos_neg[i][pos_prot1])->second).end(),(prot_domain.find(pos_neg[j][pos_prot2])->second).begin(), (prot_domain.find(pos_neg[j][pos_prot2])->second).end(),compare2.begin());
        prot1_prot4 = double(c2-compare2.begin())/((prot_domain.find(pos_neg[i][pos_prot1])->second).size()+(prot_domain.find(pos_neg[j][pos_prot2])->second).size()-double(c2-compare2.begin()));
        //prot1_prot4 = double(c2-compare2.begin())/sqrt((prot_domain.find(pos_neg[i][pos_prot1])->second).size()*(prot_domain.find(pos_neg[j][pos_prot2])->second).size());
        //prot1_prot4 = numerator2/(2-numerator2);
      }
      else {
        c2 = compare2.begin();
        prot1_prot4 = 0;
      }
      if((prot_domain.find(pos_neg[i][pos_prot2]) != prot_domain.end())&&(prot_domain.find(pos_neg[j][pos_prot1]) != prot_domain.end()))
      {
        c3 = set_intersection((prot_domain.find(pos_neg[i][pos_prot2])->second).begin(), (prot_domain.find(pos_neg[i][pos_prot2])->second).end(),(prot_domain.find(pos_neg[j][pos_prot1])->second).begin(), (prot_domain.find(pos_neg[j][pos_prot1])->second).end(),compare3.begin());
        prot2_prot3 = double(c3-compare3.begin())/((prot_domain.find(pos_neg[i][pos_prot2])->second).size()+(prot_domain.find(pos_neg[j][pos_prot1])->second).size()-double(c3-compare3.begin()));
        //prot2_prot3 = double(c3-compare3.begin())/sqrt((prot_domain.find(pos_neg[i][pos_prot2])->second).size()*(prot_domain.find(pos_neg[j][pos_prot1])->second).size());
        //prot2_prot3 = numerator3/(2-numerator3);
      }
      else {
        c3 = compare3.begin();
        prot2_prot3 = 0;
      }
      if((prot_domain.find(pos_neg[i][pos_prot2]) != prot_domain.end())&&(prot_domain.find(pos_neg[j][pos_prot2]) != prot_domain.end()))
      {
        c4 = set_intersection((prot_domain.find(pos_neg[i][pos_prot2])->second).begin(), (prot_domain.find(pos_neg[i][pos_prot2])->second).end(),(prot_domain.find(pos_neg[j][pos_prot2])->second).begin(),(prot_domain.find(pos_neg[j][pos_prot2])->second).end(),compare4.begin());
        prot2_prot4 = double(c4-compare4.begin())/((prot_domain.find(pos_neg[i][pos_prot2])->second).size()+(prot_domain.find(pos_neg[j][pos_prot2])->second).size()-double(c4-compare4.begin()));
        //prot2_prot4 = double(c4-compare4.begin())/sqrt((prot_domain.find(pos_neg[i][pos_prot2])->second).size()*(prot_domain.find(pos_neg[j][pos_prot2])->second).size());
        // prot2_prot4 = numerator4/(2-numerator4);
      }
      else {
        c4 = compare4.begin();
        prot2_prot4 = 0;
      }
#endif
      //double(c1-compare1.begin()), number of common domains, majority value are 0,1,2
#if 1
      //start:pairwise-domain (mlpk-kernel)-------------------
      //  mlpk_kernel = exp(-pow((int(c1-compare1.begin())+int(c2-compare2.begin())+int(c3-compare3.begin())+int(c4-compare4.begin()))/1,2));
      //mlpk_kernel = pow(int(c1-compare1.begin())-int(c2-compare2.begin())-int(c3-compare3.begin())+int(c4-compare4.begin()),2);
      //tppk_kernel = int(c1-compare1.begin())*int(c4-compare4.begin())+int(c2-compare2.begin())*int(c3-compare3.begin());
      // if(i == j) diag.push_back(mlpk_kernel);
      //   int a = int(c1-compare1.begin())+int(c4-compare4.begin());
      //   int b = int(c2-compare2.begin())+int(c3-compare3.begin());
      //    if(a<b) mlpk_kernel = a;
      //        else mlpk_kernel = b;
      //element += coef*mlpk_kernel;
      //element += coef*tppk_kernel;
      //end:pairwise-domain (mlpk-kernel)-------------------
#endif

      //start:domain-kernel--------------------
      //element += coef*pow(prot1_prot3+prot2_prot4-prot1_prot4-prot2_prot3,2);
      //element += coef*(prot1_prot3*prot2_prot4+prot1_prot4*prot2_prot3);
      //if((prot1_prot3+prot2_prot4)>(prot1_prot4+prot2_prot3)) domain_kernel = prot1_prot3+prot2_prot4;
      //else domain_kernel = prot1_prot4+prot2_prot3;
      double sum_prot1_prot3 = prot_domain.find(pos_neg[i][pos_prot1])->second.size()+prot_domain.find(pos_neg[j][pos_prot1])->second.size();
      double sum_prot2_prot4 = prot_domain.find(pos_neg[i][pos_prot2])->second.size()+prot_domain.find(pos_neg[j][pos_prot2])->second.size();
      double com_prot1_prot3 = double(c1-compare1.begin());
      double com_prot2_prot4 = double(c4-compare4.begin());
      double sum_prot1_prot4 = prot_domain.find(pos_neg[i][pos_prot1])->second.size()+prot_domain.find(pos_neg[j][pos_prot2])->second.size();
      double sum_prot2_prot3 = prot_domain.find(pos_neg[i][pos_prot2])->second.size()+prot_domain.find(pos_neg[j][pos_prot1])->second.size();
      double com_prot1_prot4 = double(c2-compare2.begin());
      double com_prot2_prot3 = double(c3-compare3.begin());
      double sim_prot13, sim_prot14, sim_prot23, sim_prot24;
      sim_prot13 = com_prot1_prot3/sum_prot1_prot3;
      sim_prot14 = com_prot1_prot4/sum_prot1_prot4;
      sim_prot23 = com_prot2_prot3/sum_prot2_prot3;
      sim_prot24 = com_prot2_prot4/sum_prot2_prot4;
      if(com_prot1_prot3 == 0) sim_prot13 = 0;
      if(com_prot1_prot4 == 0) sim_prot14 = 0;
      if(com_prot2_prot3 == 0) sim_prot23 = 0;
      if(com_prot2_prot4 == 0) sim_prot24 = 0;
      double sim_prot13_prot24 = sim_prot13+sim_prot24;
      double sim_prot14_prot23 = sim_prot14+sim_prot23;
      double sim_kernel;
      if(sim_prot13_prot24 > sim_prot14_prot23) sim_kernel = sim_prot13_prot24;
      else sim_kernel = sim_prot14_prot23;
     // tanimoto_sim_prot13 = double(c1-compare1.begin())/((prot_domain.find(pos_neg[i][pos_prot1])->second).size()+(prot_domain.find(pos_neg[j][pos_prot1])->second).size()-double(c1-compare1.begin())); 
      //element += coef*pow(sim_prot13_prot24-sim_prot14_prot23,2);
      element += coef*(sim_prot13*sim_prot24+sim_prot14*sim_prot23);
      //element += coef*sim_kernel;
      
      //if((double(c1-compare1.begin()) + double(c4-compare4.begin())) > (double(c2-compare2.begin()) + double(c3-compare3.begin()))) domain_kernel = double(c1-compare1.begin()) + double(c4-compare4.begin());
      //else domain_kernel = double(c2-compare2.begin()) + double(c3-compare3.begin());
      //element += coef*domain_kernel;
#if 0
      if (i == j)
      {
        domain_kernel = 1;
      }
      else{
        int max_domain_i1j1 = 0;
        if((prot_domain.find(pos_neg[i][pos_prot1])->second).size() > (prot_domain.find(pos_neg[j][pos_prot1])->second).size()){
          max_domain_i1j1 = (prot_domain.find(pos_neg[i][pos_prot1])->second).size();
        }else{
          max_domain_i1j1 = (prot_domain.find(pos_neg[j][pos_prot1])->second).size();
        }
        int max_domain_i2j2 = 0;
        if((prot_domain.find(pos_neg[i][pos_prot2])->second).size() > (prot_domain.find(pos_neg[j][pos_prot2])->second).size()){
          max_domain_i2j2 = (prot_domain.find(pos_neg[i][pos_prot2])->second).size();
        }else{
          max_domain_i2j2 = (prot_domain.find(pos_neg[j][pos_prot2])->second).size();
        }
        int max_domain_i1j2 = 0;
        if((prot_domain.find(pos_neg[i][pos_prot1])->second).size() > (prot_domain.find(pos_neg[j][pos_prot2])->second).size()){
          max_domain_i1j2 = (prot_domain.find(pos_neg[i][pos_prot1])->second).size();
        }else{
          max_domain_i1j2 = (prot_domain.find(pos_neg[j][pos_prot2])->second).size();
        }
        int max_domain_i2j1 = 0;
        if((prot_domain.find(pos_neg[i][pos_prot2])->second).size() > (prot_domain.find(pos_neg[j][pos_prot1])->second).size()){
          max_domain_i2j1 = (prot_domain.find(pos_neg[i][pos_prot2])->second).size();
        }else{
          max_domain_i2j1 = (prot_domain.find(pos_neg[j][pos_prot1])->second).size();
        }
        if((int(c1-compare1.begin()) == max_domain_i1j1)&&(int(c4-compare4.begin()) == max_domain_i2j2)){ 
          domain_kernel = 1;
        }
        else if((int(c2-compare2.begin()) == max_domain_i1j2)&&(int(c3-compare3.begin()) == max_domain_i2j1)){ 
          domain_kernel = 1;
        }
        //  else if((int(c1-compare1.begin()) == max_domain_i1j1)||(int(c2-compare2.begin()) == max_domain_i1j2)||(int(c3-compare3.begin()) == max_domain_i2j1)||(int(c4-compare4.begin()) == max_domain_i2j2)){
        //     domain_kernel = 0;
        //   }
        else domain_kernel = 0;
      }
      // element += coef*pow(prot1_prot3+prot2_prot4-prot1_prot4-prot2_prot3,2);
      element += coef*domain_kernel;
      //end:domain-kernel--------------------
#endif

#if 0
      //start:location information--------------
      vector<string> compare5(1000);
      vector<string> compare6(1000);
      vector<string> compare7(1000);
      vector<string> compare8(1000);
      vector<string>::iterator c5,c6,c7,c8;

      int location_kernel = 0;
      if((prot_location.find(pos_neg[i][pos_prot1]) != prot_location.end())&&(prot_location.find(pos_neg[j][pos_prot1]) != prot_location.end()))
      {
        c5 = set_intersection((prot_location.find(pos_neg[i][pos_prot1])->second).begin(), (prot_location.find(pos_neg[i][pos_prot1])->second).end(),(prot_location.find(pos_neg[j][pos_prot1])->second).begin(), (prot_location.find(pos_neg[j][pos_prot1])->second).end(),compare5.begin());
        //prot1_prot3 = double(c5-compare5.begin())/((prot_location.find(pos_neg[i][pos_prot1])->second).size()+(prot_location.find(pos_neg[j][pos_prot1])->second).size()-double(c5-compare5.begin()));
        //prot1_prot3 = double(c5-compare5.begin())/sqrt((prot_location.find(pos_neg[i][pos_prot1])->second).size()*(prot_location.find(pos_neg[j][pos_prot1])->second).size());
      }
      else {
        c5 = compare5.begin();
        prot1_prot3 = 0;
      }
      if((prot_location.find(pos_neg[i][pos_prot1]) != prot_location.end())&&(prot_location.find(pos_neg[j][pos_prot2]) != prot_location.end()))
      {
        c6 = set_intersection((prot_location.find(pos_neg[i][pos_prot1])->second).begin(), (prot_location.find(pos_neg[i][pos_prot1])->second).end(),(prot_location.find(pos_neg[j][pos_prot2])->second).begin(), (prot_location.find(pos_neg[j][pos_prot2])->second).end(),compare6.begin());
        //prot1_prot4 = double(c5-compare6.begin())/((prot_location.find(pos_neg[i][pos_prot1])->second).size()+(prot_location.find(pos_neg[j][pos_prot2])->second).size()-double(c6-compare6.begin()));
        //prot1_prot4 = double(c6-compare6.begin())/sqrt((prot_location.find(pos_neg[i][pos_prot1])->second).size()*(prot_location.find(pos_neg[j][pos_prot2])->second).size());
      }
      else {
        c6 = compare6.begin();
        prot1_prot4 = 0;
      }
      if((prot_location.find(pos_neg[i][pos_prot2]) != prot_location.end())&&(prot_location.find(pos_neg[j][pos_prot1]) != prot_location.end()))
      {
        c7 = set_intersection((prot_location.find(pos_neg[i][pos_prot2])->second).begin(), (prot_location.find(pos_neg[i][pos_prot2])->second).end(),(prot_location.find(pos_neg[j][pos_prot1])->second).begin(), (prot_location.find(pos_neg[j][pos_prot1])->second).end(),compare7.begin());
        //prot2_prot3 = double(c7-compare7.begin())/((prot_location.find(pos_neg[i][pos_prot2])->second).size()+(prot_location.find(pos_neg[j][pos_prot1])->second).size()-double(c7-compare7.begin()));
        //prot2_prot3 = double(c7-compare7.begin())/sqrt((prot_location.find(pos_neg[i][pos_prot2])->second).size()*(prot_location.find(pos_neg[j][pos_prot1])->second).size());
      }
      else {
        c7 = compare7.begin();
        prot2_prot3 = 0;
      }
      if((prot_location.find(pos_neg[i][pos_prot2]) != prot_location.end())&&(prot_location.find(pos_neg[j][pos_prot2]) != prot_location.end()))
      {
        c8 = set_intersection((prot_location.find(pos_neg[i][pos_prot2])->second).begin(), (prot_location.find(pos_neg[i][pos_prot2])->second).end(),(prot_location.find(pos_neg[j][pos_prot2])->second).begin(),(prot_location.find(pos_neg[j][pos_prot2])->second).end(),compare8.begin());
        //prot2_prot4 = double(c8-compare8.begin())/((prot_location.find(pos_neg[i][pos_prot2])->second).size()+(prot_location.find(pos_neg[j][pos_prot2])->second).size()-double(c8-compare8.begin()));
        //prot2_prot4 = double(c8-compare8.begin())/sqrt((prot_location.find(pos_neg[i][pos_prot2])->second).size()*(prot_location.find(pos_neg[j][pos_prot2])->second).size());
      }
      else {
        c8 = compare8.begin();
        prot2_prot4 = 0;
      }
      location_kernel = pow(int(c1-compare1.begin())-int(c2-compare2.begin())-int(c3-compare3.begin())+int(c4-compare4.begin()),2);
      //location_kernel = pow(prot1_prot3-prot1_prot4-prot2_prot3+prot2_prot4,2);
      element += coef*location_kernel;
      //end:location information--------------
#endif

#if 0
      //start:sequence information--------------
      vector<string> compare9(1000);
      vector<string> compare10(1000);
      vector<string> compare11(1000);
      vector<string> compare12(1000);
      vector<string>::iterator c9,c10,c11,c12;
      int sequence_kernel = 0;
      //k-spectrum-kernel

      vector<int> spectrum_seq1_int;
      vector<int> spectrum_seq2_int;
      map<char, int> map_aminoacid;
      int k = 3;
      int a = 0, b = 0, c = 0, d = 0;
      string seq1,seq2;
      if((prot_sequence.find(pos_neg[i][pos_prot1]) != prot_sequence.end())&&(prot_sequence.find(pos_neg[j][pos_prot1]) != prot_sequence.end())){
        seq1 = prot_sequence.find(pos_neg[i][pos_prot1])->second;
        seq2 = prot_sequence.find(pos_neg[j][pos_prot1])->second;
        //k_spectrum_kernel(seq1,seq2,k,map_aminoacid,spectrum_seq1_int,spectrum_seq2_int);
        for(int i = 0; i < spectrum_seq1_int.size(); i++){
          //ff<<spectrum_seq1_int[i]<<",";
          a += pow(spectrum_seq1_int[i]-spectrum_seq2_int[i],2);
        }
        prot1_prot3 = exp(-a);
      }
      if((prot_sequence.find(pos_neg[i][pos_prot1]) != prot_sequence.end())&&(prot_sequence.find(pos_neg[j][pos_prot2]) != prot_sequence.end())){
        seq1 = prot_sequence.find(pos_neg[i][pos_prot1])->second;
        seq2 = prot_sequence.find(pos_neg[j][pos_prot2])->second;
        //k_spectrum_kernel(seq1,seq2,k,map_aminoacid,spectrum_seq1_int,spectrum_seq2_int);
        for(int i = 0; i < spectrum_seq1_int.size(); i++){
          //ff<<spectrum_seq2_int[i]<<",";
          b += pow(spectrum_seq1_int[i]-spectrum_seq2_int[i],2);
        }
        prot1_prot4 = exp(-b);
      }
      if((prot_sequence.find(pos_neg[i][pos_prot2]) != prot_sequence.end())&&(prot_sequence.find(pos_neg[j][pos_prot1]) != prot_sequence.end())){
        seq1 = prot_sequence.find(pos_neg[i][pos_prot2])->second;
        seq2 = prot_sequence.find(pos_neg[j][pos_prot1])->second;
        //k_spectrum_kernel(seq1,seq2,k,map_aminoacid,spectrum_seq1_int,spectrum_seq2_int);
        for(int i = 0; i < spectrum_seq1_int.size(); i++){
          c += pow(spectrum_seq1_int[i]-spectrum_seq2_int[i],2);
        }
        prot2_prot3 = exp(-c);
      }
      if((prot_sequence.find(pos_neg[i][pos_prot2]) != prot_sequence.end())&&(prot_sequence.find(pos_neg[j][pos_prot2]) != prot_sequence.end())){
        seq1 = prot_sequence.find(pos_neg[i][pos_prot2])->second;
        seq2 = prot_sequence.find(pos_neg[j][pos_prot2])->second;
        //k_spectrum_kernel(seq1,seq2,k,map_aminoacid,spectrum_seq1_int,spectrum_seq2_int);
        for(int i = 0; i < spectrum_seq1_int.size(); i++){
          d += pow(spectrum_seq1_int[i]-spectrum_seq2_int[i],2);
        }
        prot2_prot4 = exp(-d);
      }
      //sequence_kernel = pow(prot1_prot3-prot1_prot4-prot2_prot3+prot2_prot4,2);
      sequence_kernel = pow(a-b-c+d,2);
      //element += coef*sequence_kernel;
      //element += sequence_kernel;
      //end:sequence information--------------
#endif

      //add phylogenetic profile
      double kernel_phy = 0;
#if 0
      if((protein_organism.find(pos_neg[i][pos_prot1]) != protein_organism.end())&&(protein_organism.find(pos_neg[j][pos_prot1]) != protein_organism.end()))
      {
        c1 = set_intersection((protein_organism.find(pos_neg[i][pos_prot1])->second).begin(), (protein_organism.find(pos_neg[i][pos_prot1])->second).end(),(protein_organism.find(pos_neg[j][pos_prot1])->second).begin(), (protein_organism.find(pos_neg[j][pos_prot1])->second).end(),compare1.begin());
        //prot1_prot3 = double(c1-compare1.begin())/((protein_organism.find(pos_neg[i][pos_prot1])->second).size()+(protein_organism.find(pos_neg[j][pos_prot1])->second).size()-double(c1-compare1.begin()));
        //prot1_prot3 = double(c1-compare1.begin())/sqrt((protein_organism.find(pos_neg[i][pos_prot1])->second).size()*(protein_organism.find(pos_neg[j][pos_prot1])->second).size());
      } else {
        c1 = compare1.begin();
        prot1_prot3 = 0;
      }
      if((protein_organism.find(pos_neg[i][pos_prot1]) != protein_organism.end())&&(protein_organism.find(pos_neg[j][pos_prot2]) != protein_organism.end()))
      {
        c2 = set_intersection((protein_organism.find(pos_neg[i][pos_prot1])->second).begin(), (protein_organism.find(pos_neg[i][pos_prot1])->second).end(),(protein_organism.find(pos_neg[j][pos_prot2])->second).begin(), (protein_organism.find(pos_neg[j][pos_prot2])->second).end(),compare2.begin());
        //prot1_prot4 = double(c2-compare2.begin())/((protein_organism.find(pos_neg[i][pos_prot1])->second).size()+(protein_organism.find(pos_neg[j][pos_prot2])->second).size()-double(c2-compare2.begin()));
        //prot1_prot4 = double(c2-compare2.begin())/sqrt((protein_organism.find(pos_neg[i][pos_prot1])->second).size()*(protein_organism.find(pos_neg[j][pos_prot2])->second).size());
      } else {
        c2 = compare2.begin();
        prot1_prot4 = 0;
      }
      if((protein_organism.find(pos_neg[i][pos_prot2]) != protein_organism.end())&&(protein_organism.find(pos_neg[j][pos_prot1]) != protein_organism.end()))
      {
        c3 = set_intersection((protein_organism.find(pos_neg[i][pos_prot2])->second).begin(), (protein_organism.find(pos_neg[i][pos_prot2])->second).end(),(protein_organism.find(pos_neg[j][pos_prot1])->second).begin(), (protein_organism.find(pos_neg[j][pos_prot1])->second).end(),compare3.begin());
        //prot2_prot3 = double(c3-compare3.begin())/((protein_organism.find(pos_neg[i][pos_prot2])->second).size()+(protein_organism.find(pos_neg[j][pos_prot1])->second).size()-double(c3-compare3.begin()));
        //prot2_prot3 = double(c3-compare3.begin())/sqrt((protein_organism.find(pos_neg[i][pos_prot2])->second).size()*(protein_organism.find(pos_neg[j][pos_prot1])->second).size());
      } else {
        c3 = compare3.begin();
        prot2_prot3 = 0;
      }
      if((protein_organism.find(pos_neg[i][pos_prot2]) != protein_organism.end())&&(protein_organism.find(pos_neg[j][pos_prot2]) != protein_organism.end()))
      {
        c4 = set_intersection((protein_organism.find(pos_neg[i][pos_prot2])->second).begin(), (protein_organism.find(pos_neg[i][pos_prot2])->second).end(),(protein_organism.find(pos_neg[j][pos_prot2])->second).begin(), (protein_organism.find(pos_neg[j][pos_prot2])->second).end(),compare4.begin());
        //prot2_prot4 = double(c4-compare4.begin())/((protein_organism.find(pos_neg[i][pos_prot2])->second).size()+(protein_organism.find(pos_neg[j][pos_prot2])->second).size()-double(c4-compare4.begin()));
        //prot2_prot4 = double(c4-compare4.begin())/sqrt((protein_organism.find(pos_neg[i][pos_prot2])->second).size()*(protein_organism.find(pos_neg[j][pos_prot2])->second).size());
      } else {
        c4 = compare4.begin();
        prot2_prot4 = 0;
      }
      //kernel_phy = pow(prot1_prot3-prot1_prot4-prot2_prot3+prot2_prot4,2);
      kernel_phy = pow(int(c1-compare1.begin())-int(c2-compare2.begin())-int(c3-compare3.begin())+int(c4-compare4.begin()),2);
      element += coef*kernel_phy;
#endif

#if 0   
      double kernel_orga;
      double a = 0;
      for(int orga = 0; orga < organism_list.size(); orga++){
        double p11 = (protein_organism.find(pos_neg[i][pos_prot1])->second)[orga];
        double p21 = (protein_organism.find(pos_neg[j][pos_prot1])->second)[orga];
        double p12 = (protein_organism.find(pos_neg[i][pos_prot2])->second)[orga];
        double p22 = (protein_organism.find(pos_neg[j][pos_prot2])->second)[orga];
        a += pow(p11-p21-p22+p12,2);
      }
      kernel_orga = exp(-a);
      element += coef*kernel_orga;

#endif

      of<<element<<" ";
      //fin<<sequence_kernel<<" ";
    }
    of<<endl;
    fin<<endl;
  }
  of.close();
  ff.close();
  fin.close();
#endif
  normalize_and_add_id(argv[feature_num+2], argv[feature_num+3], pos_neg, diag);
#endif
}
