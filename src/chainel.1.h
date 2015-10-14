#ifndef chainel_h
#define chainel_h

#include <iomanip>
#include <iostream>
#include "chainel.h"
#include <cstdlib>

using namespace std;
#define NP 5
#define COVAR

const float StepScaleFactor = 1.0;  
static  double rescale_matrix[NP*NP] = {0};  

//static unsigned long int cseed = (time(NULL)%10000); 

//#define	GRND	0.08*sqrt(-2.*log(drand48()))*cos(2*M_PI*drand48())
//#define	GRND	sqrt(-2.*log(srandom(cseed)))*cos(2*M_PI*srandom(cseed))
//#define	GRND	sqrt(-2.*log(random()))*cos(2*M_PI*random())
#define	GRND	sqrt(-2.*log(drand48()))*cos(2*M_PI*drand48())
class ChainEl{
 public:
  //start with brightest galaxy parameters
  //ChainEl(float a, float b, float c, float d, float e, float f):
  //fm_s(a), fm_i(b), fs_s(c), fs_i(d), ff_s(e), ff_i(f){};
  ChainEl(float g, float a, float b, float c, float d, float e, float f):
    cm(g), fm_s(a), fm_i(b), fs_s(c), fs_i(d), ff_s(e), ff_i(f){};
  ChainEl(){
    cm = -0.545;
    //fm_s = -0.4273;
    fm_s = -0.40;
    fm_i = 2.7458;
    fs_s = 0;
    fs_i = 2.2;
    ff_s = -0.0395;
    ff_i = 0.9128;
    SetStepsize();
  }
  float cmean(){
    return cm;
  }
  float cmean(float mr){
    return cm-0.02*(mr-1);
  }
  float csig(){
    //return 0.46;
    return 0.35;
  }
  float fmean(float mr){
    float fm = fm_s*(mr-1)+fm_i;
    if(fm<0.1) {fm = 0.1; }//cerr<<"fm prob"<<fm<<" "<<fm_s<<" "<<fm_i<<" "<<mr<<endl;}
    if(fm>7.0) {fm = 7.0;}// cerr<<"fm prob"<<fm<<" "<<fm_s<<" "<<fm_i<<" "<<mr<<endl;}
    //cout<<mr<<" "<<fm_s<<" "<<fm_i<<" "<<fm<<endl;
    return fm;
  }
  float fsig(float mr){
    float fs = fs_s*(mr-1)+fs_i;
    //if(fs<0.1) {fs = 0.1; cerr<<"fs prob"<<fs<<" "<<fs_s<<" "<<fs_i<<" "<<mr<<endl;}
    //if(fs<0.1) {fs = 0.1; cerr<<"fs prob"<<fs<<" "<<fs_s<<" "<<fs_i<<" "<<mr<<endl;}
    return fs;
  }
  float ffrac(float mr){
    float ff = ff_s*(mr-1)+ff_i;
    //    if(ff < 0.05) {ff = 0.05; cerr<<"fm prob"<<ff<<" "<<ff_s<<" "<<ff_i<<" "<<mr<<endl;}
    return ff;
  }

  int Okay(){
    //these parameters are out of bounds.  don't bother evaluating.
    int okay = -1;
    if(cm <-1.0) okay = 6;
    if(cm > 0.0) okay = 6;
    //require brighter galaxies to have lower fm
    // and require fm to be > 0 for -22.2 galaxies
    if(fm_i < 0.1) okay = 0;
    if(fm_i > 6.5) okay = 0;
    if(fm_s > 0.5) okay=1;
    if(fm_s < -1*fm_i/4+0.1) okay=11;
    //require brighter galaxies to have lower fs
    // and require fs to be > 0 for -22.2 galaxies
    if(fs_i < 0.5) okay = 2;
    if(fs_i > 4.5) okay = 2;
    if(fs_s > 0.5) okay = 3;
    if(fs_s < -1*fs_i/4) okay = 3;
    if(fs_s < -0.5) okay = 3;
    //require brighter galaxies to have lower p
    // and require p to be > 0 for -22.2 galaxies    if(ff_i < 0.1) okay = false;
    if(ff_i > 1) okay = 4;
    //    if(ff_s > (1-1*ff_i)/0.5) okay = false;
    //if(ff_s > 0.0) okay = 5;
    if(ff_s < -1*ff_i/4) okay=5;
    return okay;
  }

  void Reset(ChainEl oldel){
    //    int NP = 5;
    float buff_vec[NP];
    //float sig[NP];
    //    float sig[7]={0.17, 0.30, 0.0, 0.0, 0.085, 0.05, 0.207};
    float sig[NP]={0.207, 0.17, 0.30, 0.085, 0.05};
#ifdef COVAR
    for(int i=0;i<NP;i++)
      buff_vec[i]=GRND;//gaussian();
    int r_index;
#endif
    for(int i=0;i<NP;i++){
#ifdef COVAR               
    sig[i]=0.;
    for(int j=0;j<NP;j++){
      r_index=i*NP+j;
      sig[i]+=StepScaleFactor*rescale_matrix[r_index]*buff_vec[j];
      //      cout<<buff_vec[j]<<" "<<rescale_matrix[r_index]<<" "<<sig[i]<<endl;
      //cout<<i<<" "<<sig[i]<<endl;                       
    }
#endif
      cm =  oldel.cm + sig[0]*GRND  ;
      //this is the f_mean for an L* galaxy
      fm_s =  oldel.fm_s + sig[1]*GRND  ;
      fm_i =  oldel.fm_i + sig[2]*GRND  ;
      //this is the f_sigma for an L* galaxy
      //fs_i =  oldel.fs_i + sig[3]*GRND  ;
      //fs_s =  oldel.fs_s + sig[2]*GRND  ;
      //this is the fraction of field galaxies for L*
      ff_s =  oldel.ff_s + sig[3]*GRND  ;
      ff_i =  oldel.ff_i + sig[4]*GRND  ;
    }
  }  
    void Print()const{
      printf("%2.3f %2.3f %2.3f %2.3f %2.3f %2.3f %2.3f\t", cm, fm_s, fm_i, fs_s, fs_i, ff_s, ff_i);
    }
    
  void Write(ofstream &outfile)const{
    outfile<<setprecision(4)<<cm<<" "
	   <<fm_s<<" "<<fm_i<<" "
	   <<fs_s<<" "<<fs_i<<" "
	   <<ff_s<<" "<<ff_i<<" ";

  }


  void SetStepsize(){
    if(rescale_matrix[0] == 0){
      ifstream infile("../sdssinp/covar.dat");
      if(!infile){
	cerr<<"cannot open ../sdssinp/covar.dat";
	exit(1);
      }
      for(int i=0;i<NP*NP;i++){
	infile>>rescale_matrix[i];
      }
    }
  }






                                                                                                                        

private:
  float cm;
  float fm_s;
  float fm_i;
  float fs_s;
  float fs_i;
  float ff_s;
  float ff_i;
};
#undef GRND

#endif
