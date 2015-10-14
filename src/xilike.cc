#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "constants.h"

using namespace std;

struct  data_point {
  float x, y, err;
};

//const int np = nbin-2;

void readdata(string filename, vector <data_point> &data, int npoints){
  //  cout<<"reading"<<npoints<<endl;
  data.resize(npoints);
  ifstream infile(filename.c_str());
  if(infile.fail()){
    cerr<<"Cannot open "<<filename<<endl;
    string com = "ls -l "+filename;
    system(com.c_str());
    exit(0);
  }
  //  else cout<<"comparing to "<<filename<<endl;
  for(int i=0; i<npoints; i++){
    infile>>data[i].x>>data[i].y>>data[i].err;
  }
}

double get_chi(vector <data_point> data, vector <data_point> model){
  int nn = data.size();
  double chi2 = 0;
  ofstream out("xicomp2.out");
  //cout<<nn<<endl;
  for(int i=0; i<data.size(); i++){
    double chi2p = (data[i].y-model[i].y)*(data[i].y-model[i].y)/(data[i].err*data[i].err+model[i].err*model[i].err);
    //#define DEBUG
#ifdef DEBUG
    cout<<data[i].x<<" "<<model[i].x<<"\t";
    cout<<data[i].y<<" "<<model[i].y<<"\t";
    cout<<data[i].err<<" "<<model[i].err<<"\t"<<chi2p<<endl;
#endif
    if(data[i].x-model[i].x > 0.5){
      cout<<"[xilike] alignment problem: "<<data[i].x<<" "<<model[i].x<<" "<<chi2p<<endl;
      //      exit(1);
    }
    
    if(model[i].y > 0)
      chi2 += chi2p;
    else nn--;
  }
  //chi2 = chi2/nn;
  return chi2;
}


int xilike(string cffilename, string datfilename, double &chi2, int np){
  vector <data_point> xi_data;
  vector <data_point> xi_model;

  string dir = "../sdssinp/";
  string dfn = dir+datfilename;
  //cout<<"reading"<<dfn<<endl;
  //cout<<"reading"<<cffilename<<endl;
  readdata(dfn, xi_data, np);
  readdata(cffilename, xi_model, np);
  chi2=get_chi(xi_data, xi_model);

  return 1;
}

/*
int xilike(string nnfilename, string datfilename, double &chi2, int np){
  vector <data_point> nn_data;
  vector <data_point> _model;

  string dir = "../sdssinp/";
  string dfn = dir+datfilename;
  readdata(dfn, xi_data, np);
  readdata(cffilename, xi_model, np);
  chi2=get_chi(xi_data, xi_model);

  return 1;
}
*/
