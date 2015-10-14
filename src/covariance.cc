#include <iostream>
#include <fstream>
#include <vector>
#include "nr.h"
#include "chainel.h"
using namespace std;

void CalculateCovar(vector <ChainEl> chain, string covfile, int skipfirst){
  NRMat <double> cov1(NP,NP);  
  NRMat <double> covar(NP,NP);  
  NRVec <double> diag(NP);  
  NRVec <double> avg1(NP);  
  NRMat <double> evect(NP,NP);
  
  /*  vector <ChainEl> chain;
  ifstream infile(filename.c_str());
  NRMat <double> tmp(NP,NP);
  NRVec <double> tmp_p(7);  
  while(infile){
    int tmp;
    //expect that the chain file has two integers for each element.
    infile>>tmp>>tmp;
    for(int i=0; i<NP;i++){
      infile>>tmp_p[i];
      ChainEl chainel(tmp_p);
      chain.push_back(chainel);
    }
  }
  */

  int nsteps = chain.size();
  cout<<"calculating covariance from "<<nsteps-skipfirst<<endl;
  for(int i=0;i<NP;++i){
    avg1[i]=0;
    for(int j=0;j<NP;++j){
      cov1[i][j]=0;
      evect[i][j]=0;
    }
  }
  cout<<"cov"<<endl;
  
  for(int i=skipfirst;i<nsteps;++i){
    for(int j=0;j<NP;++j){
      avg1[j]+=chain[i][j];
      for(int k=0;k<NP;++k){
	cov1[j][k]+= chain[i][j]*chain[i][k]; 
      }
    }
  }
  
  for(int i=0;i<NP;++i)
    for(int j=0;j<NP;++j)
      covar[i][j] = cov1[i][j]/nsteps - avg1[i]*avg1[j]/nsteps/nsteps;

  NR::svdcmp(covar,diag,evect);
  //ChainEl::rotation_matrix = evect;
  //ChainEl::eigenvect = diag;
  rotation_matrix = evect;
  eigenvect = diag;
  ofstream outfile(covfile.c_str());
  for(int i=0;i<NP;++i){
    double val = 0;
    for(int j=0;j<NP;++j){
      val += evect[i][j]*sqrt(diag[j]);
    }
    outfile<<val<<"  ";
  }
  outfile<<endl;
}

