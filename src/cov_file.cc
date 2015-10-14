#include <iostream>
#include <fstream>
#include <vector>
#include "nr.h"
#include "chainel.h"
using namespace std;

void CalculateCovar(vector <ChainEl> chain, string covfile, int skipfirst){
  NRMat <double> cov1(NP,NP);  
  NRMat <double> covar(NP,NP);  
  NRMat <double> tmp(NP,NP);
  NRVec <double> tmp2(NP);
  NRVec <double> diag(NP);  
  NRVec <double> avg1(NP);  
  NRMat <double> evect(NP,NP);
  
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
    ChainEl chel = chain[i];
    if(i<10){
      chel.Print();
      cout<<endl;
    }
    for(int j=0;j<NP;++j){
      avg1[j]+=chel[j];
      for(int k=0;k<NP;++k){
	cov1[j][k]+= chel[j]*chel[k]; 
      }
    }
  }
  
  for(int i=0;i<NP;++i)
    for(int j=0;j<NP;++j)
      covar[i][j] = cov1[i][j]/nsteps - avg1[i]*avg1[j]/nsteps/nsteps;

  cout<<"covar"<<endl;
  for(int i=0;i<NP;++i){
    for(int j=0;j<NP;++j){
      cout<<covar[i][j]<<" ";
    }
    cout<<endl;
  }
  NR::svdcmp(covar,diag,evect);
  cout<<"svd"<<endl;
  for(int i=0;i<NP;++i){
    for(int j=0;j<NP;++j){
      cout<<covar[i][j]<<" ";
    }
    cout<<endl;
  }


  ofstream outfile(covfile.c_str());
  for(int i=0;i<NP;++i){
    for(int j=0;j<NP;++j){
      double val = evect[i][j]*sqrt(diag[j]);
            if((val>0)&&(val<1e-8)) val = 0;
      if((val<0)&&(val>-1e-8)) val= 0;
      outfile<<val<<"  ";
    }
    outfile<<endl;
  }
}

int main(void){
  string chainfile = "mcmcchain_bbb_hv.dat";
  //string chainfile = "mcmc_test.dat";
  string covfile = "covar_test.dat";
  ifstream infile(chainfile.c_str());
  vector <ChainEl> chain;
  while(infile){
    //double chval[NP];
    //    vector <double> chval(NP);
    NRVec <double> chval(7);  
    int tmp;
    infile>>tmp>>tmp;
    for(int i=0;i<NP; i++)
      infile>>chval[i];
    //    ChainEl chel(chval[0], chval[1], chval[2], chval[3], chval[4], chval[5], chval[6]);
    ChainEl chel(chval);
    chain.push_back(chel);
  }
  int skipfirst = 0;

  CalculateCovar(chain, covfile, skipfirst);


}
