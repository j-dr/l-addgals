#ifndef chainel_h
#define chainel_h
#include "nr.h"
#include <iomanip>
#include <iostream>
#include <cstdlib>
//#include <vector>

#define NP 7
//these are the priors
const double param_min[NP] = {-1., 0.1, -1., 0.5, -1., 0.05, -1.};
const double param_max[NP] = {0.0, 6.0, 1., 4.5, 1., 1.0, 0.1};
const bool use_param[NP] = {1, 1, 1, 1, 1, 1, 1};
const bool start_from_chain = false;
using namespace std;
const string chainstartfn = "chainstart.dat";
const string covarstartfn = "covar_aaa_hv.dat";

static NRMat <double> rotation_matrix(NP,NP);  
static NRVec <double> eigenvect(NP);
const double StepScaleFactor = 2.0/(sqrt((double) NP));  
//static  double rescale_matrix[NP*NP] = {0};  

static int cseed = 100;//(time(NULL)%10000); 
class ChainEl{
 public:
  //start with brightest galaxy parameters
  //  ChainEl(double & chval){
  //  ChainEl(vector <double> chval){
  ChainEl(float g, float a, float b, float c, float d, float e, float ff){
    chainval[0] = g;
  }


  ChainEl(NRVec <double> chval){
    for(int i = 0; i<NP; i++)
      chainval[i] = chval[i];
    
  }
  ChainEl(){
    cseed = 1000;
    ifstream infile(chainstartfn.c_str());
    for(int i = 0; i<NP; i++){
      infile>>chainval[i];
    }
    cout<<endl;
    if(start_from_chain)
      SetStepSize(covarstartfn.c_str());
    else SetInitialSteps();
    PrintRotation();
  }
  //  double operator[](int index) const;
  double operator[](int index) const{
    double return_me=0;
    //    cout<<index<<"index"<<endl;
    assert( (index >= 0) && (index <NP) ); // bounds checking
    return_me = chainval[index];
    return return_me;
  }
  
  double cmean(){
    double cm  = chainval[0];
    return cm;
  }
  double csig(){
    return 0.46;
  }
  double fmean(double mr){
    double fm_i = chainval[1];
    double fm_s = chainval[2];
    double fm = fm_s*(mr-1)+fm_i;
    if(fm<0.1) {fm = 0.1; }
    if(fm>7.0) {fm = 7.0;}
    return fm;
  }
  double fsig(double mr){
    double fs_i = chainval[3];
    double fs_s = chainval[4];
    double fs = fs_s*(mr-1)+fs_i;
    return fs;
  }
  double ffrac(double mr){
    double ff_i = chainval[5];
    double ff_s = chainval[6];
    double ff = ff_s*(mr-1)+ff_i;
    if(ff < 0.05) ff = 0.05;
    if(ff > 1 ) ff = 1.0;
    return ff;
  }
  int Okay()const;
  void Write(ofstream &outfile)const;
  void Reset(ChainEl oldel);
  void Print()const;
  void PrintRotation()const;
  void SetStepSize(string covfile);

  void SetInitialSteps(){
    for(int i=0;i<NP;i++){
      rotation_matrix[i][i] = 1;
      eigenvect[i] = 0.02*(param_max-param_min);
    }

  }
  
 private:
  double chainval[NP];
  /* double cm;
     double fm_i; 
     double fm_s;
     double fs_i;
     double fs_s;
     double ff_i;
     double ff_s;
  */
};

#endif
