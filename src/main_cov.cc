#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cstdlib>
#include <cassert>
#include <unistd.h>     //for sleep()                                                   
#define NP 7
#include "nr.h"
#include "cosmo.h"
#include "constants.h"


using namespace std;

class ChainEl{
public:
  //  ChainEl(double a, double b, double c, double d):aa(a), bb(b), cc(c), dd(d){};
  //  ChainEl(vectocor <double> a){
  ChainEl(NRVec <double> a){
    for(int i=0;i<NP;i++)
      chainval[i] = a[i];
  }
  double operator[](int index) const{
    double return_me=0;
    //    cout<<index<<"index"<<endl;
    assert( (index >= 0) && (index <NP) ); // bounds checking
    return_me = chainval[index];
    return return_me;
  }
private:
  double chainval[NP];
};

//double rescale_matrix[NP*NP] = {0};  
NRMat <double> cov1(NP,NP);  
NRMat <double> covar(NP,NP);  
NRMat <double> tmp(NP,NP);
NRMat <double> tmp1(NP,NP);
NRVec <double> tmp2(NP);
NRVec <double> diag(NP);  
NRVec <double> avg1(NP);  
NRMat <double> evect(NP,NP);
int nrot;

int main(void){
  int seed = 1000;
  
  for(int i=0; i<10000; i++){
    float a = NR::gasdev(seed);
    cout<<a<<endl;
  }		
  return 0;
  
  ifstream infile("chain.dat");

  std::vector <ChainEl> chain;
  
  while(infile){
    //for(int i=0; i<NSTEP;i++){
    //    vector <double> p(NP);
    for(int i=0;i<NP; i++)
      infile>>tmp2[i];
    ChainEl chainel(tmp2);
    chain.push_back(chainel);
  }
  int NSTEP = chain.size();
    cout<<NSTEP<<endl;
    for(int j=0;j<NP;++j){
    avg1[j]=0;
    for(int k=0;k<NP;++k){
      cov1[j][k]=0;
    evect[j][k]=0;
    }
  }
  cout<<"cov"<<endl;

  for(int i=0;i<NSTEP;++i){
    for(int j=0;j<NP;++j){
      avg1[j]+=chain[i][j];
      for(int k=0;k<NP;++k){
	//cout<<i<<" "<<j<<" "<<k<<endl;
	cov1[j][k]+= chain[i][j]*chain[i][k]; 
	//	cout<<avg1[j]<<" "<<cov1[j][k]<<endl;
      }
    }
  }

  for(int i=0;i<NP;++i){
    for(int j=0;j<NP;++j){
      covar[i][j] = cov1[i][j]/NSTEP - avg1[i]*avg1[j]/NSTEP/NSTEP;
      //   tmp[i][j] = cov1[i][j]/NSTEP - avg1[i]*avg1[j]/NSTEP/NSTEP;
      cout<<covar[i][j]<<' ';
    }
    cout<<endl;
  }
  /*
  ifstream infile2("jeremy2.dat");
  
  for(int i=0;i<NP;++i){
    for(int j=0;j<NP;++j){
      infile2 >>covar[i][j];
    }
  }

  cout<<"jeremy cov"<<endl;
  for(int i=0;i<NP;++i){
    for(int j=0;j<NP;++j){
      cout<<covar[i][j]<<' ';
    }
    cout<<endl;
  }
  */

  NR::jacobi(covar,diag,evect,nrot);
  cout<<"jacobi"<<endl;
  for(int i=0;i<NP;++i){
    for(int j=0;j<NP;++j){
       cout<<evect[i][j]<<' ';
    }
    cout<<endl;
  }

  cout<<"step 1"<<endl;
  for(int i=0;i<NP;++i){
    tmp2[i] = 0;
    for(int j=0;j<NP;++j){
	tmp2[i] += evect[i][j]*sqrt(diag[j]);
    }
    cout<<tmp2[i]<<endl;
  }


/*
  for(int i=0;i<NP;++i){
    cout<<diag[i]<<' ';
    cout<<endl;
  }


  NR::gaussj(evect,tmp1);

  //double fac=1;//2.4/(sqrt(NP));

  cout<<"gaussj"<<endl;
  for(int i=0;i<NP;++i){
    for(int j=0;j<NP;++j){
      cout<<evect[i][j]<<' ';
    }
    cout<<endl;
  }
*/
  NR::svdcmp(covar,diag,evect);
  cout<<"svd  -- u "<<endl;
  for(int i=0;i<NP;++i){
    for(int j=0;j<NP;++j){
       cout<<evect[i][j]<<' ';
    }
    cout<<endl;
  }

  cout<<"step 2"<<endl;
  for(int i=0;i<NP;++i){
    tmp2[i] = 0;
    for(int j=0;j<NP;++j){
	tmp2[i] += evect[i][j]*sqrt(diag[j]);
    }
    cout<<tmp2[i]<<endl;
  }

  return 0;
  /*

  NR::svdcmp(tmp,diag,evect);

  cout<<"svd  -- u "<<endl;
  for(int i=0;i<NP;++i){
    for(int j=0;j<NP;++j){
       cout<<evect[i][j]<<' ';
    }
    cout<<endl;
  }


  cout<<"svd -- v"<<endl;
  for(int i=0;i<NP;++i){
    for(int j=0;j<NP;++j){
       cout<<evect[i][j]<<' ';
    }
    cout<<endl;
  }


  cout<<"svd -- vT"<<endl;
  for(int i=0;i<NP;++i){
    for(int j=0;j<NP;++j){
       cout<<evect[j][i]<<' ';
    }
    cout<<endl;
  }

  cout<<"svd = [ortho_matrix]T"<<endl;
  for(int i=0;i<NP;++i){
    for(int j=0;j<NP;++j){
      cout<<evect[i][j]<<' ';
    }
    cout<<endl;
  }



  for(int i=0;i<NP;++i){
    for(int j=0;j<NP;++j){
      //      cout<<evect[j][i]*sqrt(diag[i])*evect[i][j]<<' ';
      //cout<<sqrt(diag[j])*evect[i][j]<<' ';
      tmp[i][j] = sqrt(diag[j])*evect[i][j];
    }
    //    cout<<endl;
  }

  cout<<"sqrt(diag)*eigenvectors"<<endl;
  for(int i=0;i<NP;++i){
    for(int j=0;j<NP;++j){
      tmp[i][j] = sqrt(diag[j])*evect[i][j];
      //      for(int k=0;k<NP;++k){
      //tmp[j][i];
      cout<<tmp[i][j]<<" ";
    }
    cout<<endl;
  }
  */
  /*  cout<<"svd"<<endl;
  for(int i=0;i<NP;++i){
    for(int j=0;j<NP;++j){
      //      cout<<evect[j][i]*sqrt(diag[i])*evect[i][j]<<' ';
      cout<<evect[i][j]*tmp[i][j]<<" ";
    }
    cout<<endl;
  }

  */

  cout<<"matrix 1"<<endl;
  for(int i=0;i<NP;++i){
    for(int j=0;j<NP;++j){
      tmp1[i][j] = 0;
      for(int k=0;k<NP;++k){
	//	cout<<i<<" "<<j<<" "<<k<<" "<<evect[k][j]<<" "<<tmp[i][k]<<endl;
	//	tmp1[i][j] += evect[k][j]*tmp[i][k];
	//	tmp1[i][j] += evect[j][k]*tmp[k][i];
	//tmp1[i][j] += evect[i][k]*tmp[k][j];
	tmp1[i][j] += evect[j][k]*tmp[k][i];
      }
      cout<<tmp1[i][j]<<" ";
    }
    cout<<endl;
  }

  /*
  cout<<"matrix 2"<<endl;
  for(int i=0;i<NP;++i){
    for(int j=0;j<NP;++j){
      tmp1[i][j] = 0;
      for(int k==0;k<NP;++k)
	tmp1[k][j] += evect[k]*tmp[i];
      cout<<tmp1[i][j]<<" ";
    }
    cout<<endl;
  }
  */
  
      return 0;
}



