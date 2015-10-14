/*

calculate_wp function takes input file of xi(r) data
and calculates wp(r) in specified bins (not required
to be the same as xi(r).  extrapolates a power-law 
beyond xi(rmax) and below xi(rmin).

Risa Wechsler Jan 2005

*/

#include "utils.h"
#include "integrator.h"
#define sqr(x) ((x)*(x))
/*#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>



using namespace std;

struct data_point {
  float x, y, err;
};


void read_data(string filename, data_point data[], int npoints);
//void read_cdata(string filename, vector < vector <float> > &data, int npoints);
void	spline(float x[], float y[] ,int n, float yp1, float ypn, float y2[]);
void	splint(float xa[], float ya[], float y2a[], int n, float x, float *y);
void	spldint(float xa[], float ya[], float y2a[], int n, float x, float *y);
double	gammafn2(float a, float xmin);
double	gammln(float xx);
double	gammafn1(float z);
*/




//#define NBINS 11
#define NBINS2 20

static float rrr[NBINS2];
static float xxx[NBINS2];
static float splxi[NBINS2];

static double xirp(double y, double* params){
  float xi;
  float rp = params[0];
  float rf = sqrt(sqr(rp)+sqr(y)) ;
  float logam = params[1];
  float higam = params[2];
  float r0 = params[3];
  int nbins = (int) (params[4]);
  //cout<<"*"<<nbins<<endl;
  float logxir = 0;
  if(log10(rf)<rrr[0]){
    xi  = pow(rf/r0,logam);
  }
  else if((log10(rf)<rrr[NBINS2-1])){
    splint(rrr, xxx, splxi, NBINS2, log10(rf), &logxir);
    xi = 2.*pow((float) (10.),logxir);
  }
  else{
    xi  = pow(rf/r0,higam);
  }
  return xi;
}

float wp(float rp, float ymax, float logam, float higam, float r0, int nb){
  integrator xirp_int(&xirp, 5, 1.e-5, .01);
  return xirp_int.Integrate(0, ymax, rp, logam, higam, r0, nb);
}

int calculate_wp(string fname, string outfname, int nbins)
{
  float ymax=80.0 ;
  //  float dy=0.005;
  
  data_point cfdata[NBINS2];
  //read_data(fname, cfdata, nbins);
  read_data(fname, cfdata, NBINS2);
  int j = 0;
  for(int i=0;i<NBINS2;i++){
    //cout<<"[wp]"<<i<<" "<<j<<" "<<cfdata[i].x<<" "<<cfdata[i].y<<" "<<nbins<<endl;
    if(cfdata[i].y>-0.99){ //don't save bad values (xi(r) = -1)
      rrr[j] = log10(cfdata[i].x);
      xxx[j] = log10(cfdata[i].y);
      cout<<rrr[j]<<" "<<xxx[i]<<endl;
      j++;
    }
  }
  cout<<"done"<<endl;
  nbins = j; //reset to ignore bad values

  float y1 = (xxx[0]-xxx[1])/(rrr[0]-rrr[1]);
  float yn = (xxx[NBINS2-2]-xxx[NBINS2-1])/(rrr[NBINS2-2]-rrr[NBINS2-1]);
  spline(rrr, xxx, NBINS2, y1, yn, splxi);

  //float rmin = 0.257;
  float dlogr = 0.2;
  //float rmax = 16.2181;

  // note that this is a somewhat crude approximation, but the
  // accuracy is enough not to affect the result within the scales we 
  // are intersted in 
  float logam = (xxx[0]-xxx[3])/(rrr[0]-rrr[3]);
  float higam = (xxx[NBINS2-4]-xxx[NBINS2-1])/(rrr[NBINS2-4]-rrr[NBINS2-1]);
  float gam = (xxx[1]-xxx[NBINS2-2])/(rrr[1]-rrr[NBINS2-2]);
  float r0 = pow(10, -1*(xxx[NBINS2-3]-gam*rrr[NBINS2-3])/gam);
  cout<<logam<<" "<<higam<<" "<<gam<<" "<<r0<<endl;
  ofstream outfile(outfname.c_str());
  for(int j=0; j<nbins;j++){
    float logrmin = -0.79;//log10(
    //  for(float logr = logrmin; logr
    //  float rp = pow(10, log10(rmin)+j*dlogr);
    float rp = pow(10, logrmin+j*dlogr);
    //    cout<<j<<" "<<rp<<endl;
    float wpgg = wp(rp, ymax, logam, -1.64, r0, nbins);
    outfile<<rp<<" "<<wpgg<<" "<<0<<endl;
  }
  return 0 ;
}
