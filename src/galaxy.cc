#include <cmath>
#include <iostream>
#include "nr.h"
#include "myrand.h"
#include "fivetuple.h"
#include "constants.h"
#include "galaxy.h"
#include "cosmo.h" 
#include "ReadParameters.h"
#include "owens.hpp"
#ifdef LF_FROM_DATA
#include <vector>
#endif

using namespace std;

void read_growthfcn2(vector <float> zarr, vector <float> darr){
  std::string filename = "/afs/slac.stanford.edu/u/ki/mbusha/projects/modules/idl/cosmology/growthfcn_tab.ascii";
  ifstream file(filename.c_str());
  if (file.fail()) {
    std::cerr<<"error: cannot open "<<filename<<endl;
    exit(1);
  }
  else std::cout<<"reading "<<filename<<endl;
  int entries=2000;
  zarr.resize(entries);
  darr.resize(entries);
  for (int ii=0; ii<entries; ii++) {
    file>>zarr[ii]>>darr[ii];
  }
  file.close();
}

float FiveTuple::LocalDens() const{

  // returns a point selected from the
  // sum of a gaussian and a lognormal
  // with parameters specifed by the ftuple
  float min = -2.30259; //log(0.1)
  float max = 2.14007;  //log(8.5)
  float d8;
  if(randbool(p)){ //field galaxy
    float f8=-1;
    if((fm<0.1)||(fm>8.5)) cout<<"extreme f8"<<fm<<" "<<fs<<endl;
    while((f8<0.1)||(f8>8.5)) 
      f8 = fm+sqrt(fs)*NR::gasdev(seed);
    d8 = f8;
  }
  else{ //central galaxy
    float c8 = min;
    if(cm<min) cout<<"extreme c8"<<cm<<" "<<cs<<endl;
    while((c8<=min)||(c8>max))
      c8 = cm+sqrt(cs)*NR::gasdev(seed);
    d8 = exp(c8);
  }
  return d8;
}

//NBIN is the total number of rdel bins for our magnitude calculation
//#define NBIN 85000
#define NBIN 8500
struct den_ent{
  float prob[NBIN];
  float r[NBIN];
};

float LocalDens(FiveTuple fp, FiveTuple fpp, int weight1, int weight2)
{
  //generate a table of P(y) = x
  cout<<"in LocalDens..."<<endl;
  float x[NBIN];
  float y[NBIN];
  int nbin = NBIN;
  float d8_max = 8.5;
  float y_del = d8_max/nbin;
  for(int i=0;i<nbin;i++){
    y[i] = (i+1)*y_del;
    float p1 = 0.5*(1. - fp[4])*(1+erf((log(y[i])-fp[0])/(fp[1]*sqrt(2.0))));
    float p2 = 0.5*fp[4]*(1+erf((y[i]-fp[2])/(fp[3]*sqrt(2.0))));
    float p3 = 0.5*(1. - fpp[4])*(1+erf((log(y[i])-fpp[0])/(fpp[1]*sqrt(2.0))));
    float p4 = 0.5*fpp[4]*(1+erf((y[i]-fpp[2])/(fpp[3]*sqrt(2.0))));
#ifdef SIXPARAMS
    cout<<t((y[i]-fp[2])/fp[3],fp[5])<<endl;
    p2 -= fpp[4]*2.0*t((y[i]-fp[2])/fp[3],fp[5]);
    p4 -= fpp[4]*2.0*t((y[i]-fpp[2])/fpp[3],fpp[5]);
#endif
    x[i] = weight1*(p1+p2) - weight2*(p3+p4);
  }
  //normalize to make sure it equals one -- easy since we're working with the integral
  float norm = x[nbin-1];
  for(int i=0;i<nbin;i++)
    x[i] /= norm;
  float ranu = drand48();
  int ind = 0;
  for(ind=0;ind<nbin;ind++){
    if(ranu < x[ind])
      break;
  }
  if (ind == 0)
    ind++;
  if (ind == nbin)
    //nbin++;
    ind--;
  float dx = ranu - x[ind-1];
  float slope = (y[ind]-y[ind-1])/(x[ind]-x[ind-1]);
  float d8 = y[ind-1]+dx*slope;
  if (d8 > d8_max)
    d8 = d8_max;
  return d8;
}

float LocalDens(den_ent pdf)
{
  float d8_max = 8.5;
  //float d8_min = 0.01;
  float d8_min = 1.1*minrnn;
  int nbin = NBIN;
  float d8 = -1.0;
  while(d8 < d8_min || d8 > d8_max){
    float ranu = drand48();
    int ind = 0;
    for(ind=0;ind<nbin;ind++){
      if(ranu < pdf.prob[ind])
	break;
    }
    if (ind == 0)
      ind++;
    if (ind == nbin)
      //nbin++;
      ind--;
    float dx = ranu - pdf.prob[ind-1];
    float slope = (pdf.r[ind]-pdf.r[ind-1])/(pdf.prob[ind]-pdf.prob[ind-1]);
    d8 = pdf.r[ind-1]+dx*slope;
  }
  if(!(d8 < d8_max)) d8 = d8_max;
  if(!(d8 > d8_min)) d8 = d8_min;
  return d8;
}

float SelectGalaxyZ()
{
  //selects a random radius for galaxy s.t. prob \propto r^3
  float rn = drand48();
  float rm = cosmo.RofZ(ZREDMAX);
  float rmin = cosmo.RofZ(ZREDMIN);
  //rn *= rm*rm*rm;
  rn = (rm*rm*rm - rmin*rmin*rmin)*rn + rmin*rmin*rmin;
  rn = pow(rn, 1.0/3.0);
  float z = cosmo.ZofR(rn);
  return z;
}

float SelectGalaxyZ(float rmax)
{
  //selects a random radius for galaxy s.t. prob \propto r^3
  float rn = drand48();
  float rm = rmax;
  //rn *= rm*rm*rm*rm;
  //rn = pow(rn, 0.25);
  rn *= rm*rm*rm;
  rn = pow(rn, 1.0/3.0);
  float z = cosmo.ZofR(rn);
  return z;
}

#ifdef LF_FROM_DATA
vector <float> gmagbin;
vector <float> gdensity;

void read_lf_data(void){
  string lf_file = "LF.dat";

  //float rho = sim.Np()*sim.ParticleMass()/volume;
  //float rhoback = 3.*100*100/(8*!PI*4.301e-9)*sim.OmegaM();
  //float phi_rescale = rho/rhoback;
  //cout<<"Rescaling factor for LF normalization: "<<phi_rescale<<endl;
  //phi_rescale = 1;

  float gmagbin1, gdensity1;
  ifstream file(lf_file.c_str());
  cout<<"Reading LF file "<<lf_file<<endl;
  if (file.fail()){
    cerr<<"[galaxy.cc] ERROR:  Can't open lf file: `"<<lf_file<<"' for reading."<<endl;
    exit(1);
  }
  while(file){
    file>>gmagbin1>>gdensity1;
    gmagbin.push_back(gmagbin1);
    gdensity.push_back(gdensity1i);
  }
  file.close();
}
#endif

double MLF(double M, double* dummy){
  // return pow(10,-0.4*(M-Mstar))*0.4*log(10.0)*phistar*pow(10,-0.4*(M-Mstar)*(alpha+1))
  // *exp(-1*pow(10,-0.4*(M-Mstar)));
 return (M-Mstar)*0.4*log(10.0)*phistar*pow(10.,-0.4*(M-Mstar)*(alpha+1))
   *exp(-1*pow(10.,-0.4*(M-Mstar)));

}

double LF(double M, double* dummy){
  //SDSS luminosity function in the I band
  //from Blanton et al 02
  //  double Mstar = -20.82;
  // double phistar = 0.0147;
  // double alpha = -1.00;

#ifdef LF_FROM_DATA
 //We're using a tabulated LF read in from a file
  int i = 0;
  while(M > gmagbin[i] && i < gmagbin.size())
    i++;
  if (i == 0)
    i++;
  float diff = M - gmagbin[i-1];
  float LF,slope;
  if (gdensity[i] == 0 || gdensity[i-1] == i){
    slope = (gdensity[i]-gdensity[i-1])/(gmagbin[i]-gmagbin[i-1]);
    LF = gdensity[i-1] + diff*slope;
    if (LF <= 0.)
      LF = 0.;
  } else {
    slope = log10(gdensity[i]/gdensity[i-1])/(gmagbin[i]-gmagbin[i-1]);
    LF = pow(10., log10(gdensity[i-1]) + diff*slope);
  }
  return LF;
#else
  //Just use a schechter LF
  return 0.4*log(10.0)*phistar*pow(10.,-0.4*(M-Mstar)*(alpha+1))
    *exp(-1*pow(10.,-0.4*(M-Mstar)));
#endif
}

float ChooseMag(){
  unsigned int sample_size = 1000;
  vector <double> mags;
  GetMags(sample_size,mags);
  //  int s = randint(&seed, sample_size);
  int s = randint(sample_size);
  return mags[s];
}

void GetMags(unsigned int n, vector <double> &mags){
  integrator LF_integrator(&LF, 0, 1.0E-05, 0.01);
  mags.reserve(n);
  double volume = cube(3000/16.*use_cells);//n/0.0370370;
  //cout<<"volume: "<<volume<<" "<<cube(3000/16.)<<endl;
  for(double mag=-23.5; mag <=-8; mag+=0.0001){
      double num = LF_integrator.Integrate(-30, mag)*volume;
      int diff = ((int) floor(num))-mags.size();
      if (num>n) diff = n-mags.size();
      for(int i=0;i<diff;i++)
	mags.push_back(mag);
      if (mags.size()>=n) break;
  }  
}

void GetMags(double vol, vector <double> &mags){
  integrator LF_integrator(&LF, 0, 1.0E-05, 0.01);
  int n = (int) (LumNumberDensity(Magmin)*vol);
  mags.reserve(n);
  double volume = vol;
  cout<<volume<<" "<<cube(3000/16.*use_cells)<<endl;
  for(double mag=-23.5; mag <=-8; mag+=0.0001){
      double num = LF_integrator.Integrate(-30, mag)*volume;
      int diff = ((int) floor(num))-mags.size();
      if (num>n) diff = n-mags.size();
      for(int i=0;i<diff;i++)
	mags.push_back(mag);
      if (mags.size()>=n) break;
  }  
}

void ReadLFFile(void){
  string filename = "../sdssinp/lumnumberdens.dat";
  /*
#ifndef LF_FROM_DATA
  ifstream file(filename.c_str());
  if (file.fail()) {
    cerr<<"error: cannot open "<<filename<<endl;
    exit(1);
  }
  else cout<<"reading "<<filename<<endl;
  //int entries=3801;
  int entries=3201;
  magnitude.resize(entries);
  lumnumdens.resize(entries);
  for (int ii=0; ii<entries; ii++) {
    file>>magnitude[ii]>>lumnumdens[ii];
  }
  file.close();
#else
  */
  int entries=3201;
  double* dummy;
  magnitude.resize(entries);
  lumnumdens.resize(entries);
  for (int ii=0; ii<entries; ii++) {
    //lumnumdens[ii] = -24+ii*0.005;
    //magnitude[ii] = LF(lumnumdens[ii], dummy);
    magnitude[ii] = -24+ii*0.005;
    lumnumdens[ii] = LF(magnitude[ii], dummy);
  }
  //#endif
}


double LumNumberDensityInterp(double M){
  //if((M<-25)||(M>-5)){
  if((M<-30)||(M>5)){
    cerr<<"Magnitude out of bounds in LumNumberDensity function: "<<M<<endl;
    exit(1);
  }
  vector<double>::iterator pos;
  pos=upper_bound(magnitude.begin(), magnitude.end(), M);

  int i = distance(magnitude.begin(),pos);
  if (i >= magnitude.size() - 2) 
    i = magnitude.size() - 3;
  if (i <= 0) i = 1;
  double delta_x = magnitude[i] - magnitude[i-1];
  double dx      = M        - magnitude[i-1];
  double delta_y = lumnumdens[i] - lumnumdens[i-1];
  double y = dx/delta_x*delta_y+lumnumdens[i-1];

//cout<<" closest mag index, value: "<<i<<" "<<magnitude[i]<<endl;
//cout<<" delta_x, dx: "<<delta_x<<" "<<dx<<endl;
//cout<<" delta_y, y: "<<delta_y<<" "<<y<<endl;

  /*
  int i = distance(magnitude.begin(),pos);
  i--;
  double delta_x = magnitude[i] - magnitude[i-1];
  double dx      = M        - magnitude[i];
  double delta_y = lumnumdens[i] - lumnumdens[i-1];
  double y = dx/delta_x*delta_y+lumnumdens[i-1];
//cout<<"in Interp:  M="<<M<<", i="<<" "<<i<<", delta_x="<<delta_x<<", dx="<<dx<<", delta_y="<<delta_y<<", lumnumdnes[i-1]="<<lumnumdens[i-1]<<", y="<<y<<endl;
  */

  return y;
}



void ReadDMFile(void){
  string filename = "../sdssinp/zdistmod.dat";
  ifstream file(filename.c_str());
  if (file.fail()) {
    cerr<<"error: cannot open "<<filename<<endl;
    exit(1);
  }
  else cout<<"reading "<<filename<<endl;
  int entries=749;
  zdistmod_z.resize(entries);
  zdistmod_dm.resize(entries);
  for (int ii=0; ii<entries; ii++) {
    file>>zdistmod_z[ii]>>zdistmod_dm[ii];
  }
  cout<<"reading "<<filename<<" "<<endl;
  file.close();
}


double ZdistmodInterp(double z){
  if((z<0)||(z>1.5)){
    cerr<<"Redshift out of bounds in Zdistmod function: "<<z<<endl;
    exit(1);
  }
  vector<double>::iterator pos;
  pos=upper_bound(zdistmod_z.begin(), zdistmod_z.end(), z);
  int i = distance(zdistmod_z.begin(),pos);
  double delta_x = zdistmod_z[i] - zdistmod_z[i-1];
  double dx      = z        - zdistmod_z[i-1];
  double delta_y = zdistmod_dm[i] - zdistmod_dm[i-1];
  double y = dx/delta_x*delta_y+zdistmod_dm[i-1];

  return y;
}

double LumNumberDensity(double M){
  integrator LF_integrator(&LF, 0, 1.0E-05, 0.01);
  return LF_integrator.Integrate(-30, M);
}

double LumNumberDensity(double M1, double M2){
  integrator LF_integrator(&LF, 0, 1.0E-05, 0.01);
  return LF_integrator.Integrate(M1, M2);
}


double AverageMagnitude(double M1, double M2){
  integrator LF_integrator(&LF, 0, 1.0E-05, 0.01);
  float tot = LF_integrator.Integrate(M1, M2);
  for(double mag=M1; mag<=M2;mag+=0.001){
   float thistot = LF_integrator.Integrate(M1, mag);
   if (thistot>0.5*tot){
     return mag;
     break;
   }
  }
  return 0;
}

/*double NdensLum(double ndens){
  double mag =  NdensMagnitude(ndens);
  double loverlstar = pow(10,-0.4*(mag-Mstar));
  return loverlstar*Lstar;
  }*/

//Returns the magnitude of a specified number density
double NdensMagnitude(double ndens){
  double mmin = -30.;
  double mmax = -8.;
  double diff = 1.;
  double nndensmid;
  double mmid;
  //cout<<ndens<<" "<<LumNumberDensity(mmin)<<" "<<LumNumberDensity(mmax)<<" ";
  //cout<<magnitude.size()<<" "<<magnitude[0]<<" "<<magnitude[magnitude.size()-1]<<" "<<lumnumdens[0]<<" "<<lumnumdens[magnitude.size()-1]<<endl;cout<<magnitude[0]<<" "<<magnitude[magnitude.size()-1]<<" "<<lumnumdens[0]<<" "<<lumnumdens[magnitude.size()-1]<<endl;
  while(diff>0.0001){
    mmid=(mmin+mmax)/2.;
    //    cout<<"range"<<ndens<<" "<<mmin<<" "<<mmax<<endl;
    double nndens1 = LumNumberDensityInterp(mmin);
    double nndens2 = LumNumberDensityInterp(mmax);
    nndensmid = LumNumberDensityInterp((mmin+mmax)/2);
    //if((ndens<nndens1)||(ndens>nndens2)){
    //cout<<mmin<<" "<<mmid<<" "<<mmax<<endl;
    //cout<<"bad val: "<<nndens1<<" "<<ndens<<" "<<nndens2<<endl;
    //}
    if(ndens<nndens1){
      cout<<mmin<<" "<<mmid<<" "<<mmax<<endl;
      cout<<"bad val: "<<nndens1<<" "<<ndens<<" "<<nndens2<<endl;
    }
    if(ndens>nndens2){
      cout<<mmin<<" "<<mmid<<" "<<mmax<<endl;
      cout<<"bad val: "<<nndens1<<" "<<ndens<<" "<<nndens2<<endl;
    }
    //assert((ndens>nndens1)&&(ndens<nndens2));
    if(!(ndens>=nndens1)){
      cout<<"Going to fail assert for nndens1!"<<endl;
      cout<<"bad val: "<<nndens1<<" "<<ndens<<" "<<nndens2<<endl;
      cout<<"Mags: "<<mmin<<" "<<mmax<<endl;
    }
    assert(ndens>=nndens1);
    if(!(ndens<=nndens2)){
      cout<<"Going to fail assert for nndens2!"<<endl;
      cout<<"bad val: "<<nndens1<<" "<<ndens<<" "<<nndens2<<endl;
      cout<<"Mags: "<<mmin<<" "<<mmax<<endl;
    }
    assert(ndens<=nndens2);
    diff = mmax-mmin;
    if(diff<0.0001) break;
    if(ndens>nndensmid) mmin=mmid;
    else mmax=mmid;
  }

  return mmid;
}

/*
  Short routine to readin the denspdf parameters.  Assumes that the 
  input is an ascii with 6 columns, denoting magnitude, cm, cs, fm, fs, p
*/
void ReadPDFs(vector <FiveTuple> &denspdf, vector <float> &mags){

  ifstream file(denspdffile.c_str());
  if (file.fail()) {
    cerr<<"error: cannot open "<<denspdffile<<endl;
    exit(1);
  }
  else
    cerr<<"reading "<<denspdffile<<endl;

  float mag, cm, cs, fm, fs, p, alpha;
  int nparams = 4;

  cout<<"Reading file '"<<denspdffile<<"' with magniutude + 5 parameter PDF"<<endl;
  nparams = 5;
  while(file){
#ifdef SIXPARAMS
    file>>mag>>cm>>cs>>fm>>fs>>p>>alpha;
    FiveTuple fTup(cm,cs,fm,fs,p,alpha);
#else
    file>>mag>>cm>>cs>>fm>>fs>>p;
    FiveTuple fTup(cm,cs,fm,fs,p);
#endif
#ifdef DEBUG
    fTup.Print();
#endif
    mags.push_back(mag);
    denspdf.push_back(fTup);
  }   

  cout<<"Read "<<denspdf.size()<<" lines from denspdf"<<endl<<endl;
}

/*
  Subroutine that returns the Probability for a galaxy with magnitude M to have local density rdel
  by converting from our magnitude limited denspdf to one in a magnitude bin.  Requires information 
  on the global luminosity function
*/
void define_prob(vector <den_ent> &pdf, vector <int> &weight, vector <FiveTuple> &denspdf, vector <float> &dmagbins, float vol)
{
  //fill "weight" with the total number of galaxies brighter than dmagbins[i]
  for(int i=0;i<dmagbins.size();i++){
    int tweight = (int) (LumNumberDensity(dmagbins[i])*vol);
    weight.push_back(tweight);
  }  

  //Note that NBIN defines the total number of probability bins we'll be using in r
  int nbin = NBIN;
  float rdel = 8.5/nbin;

  //Reminder:  denspdf.size is the total number of magnitude bins
  cout<<"calculating the probabilities..."<<endl;
  for(int id=0;id<denspdf.size()-1;id++){
    //The difference in weight1 and weight2 gives the # of galaxies in this mag bin
    float weight1 = weight[id];
    float weight2 = weight[id+1];

    //The CUMULATIVE density PDFs both above and below our current magnitude bin
    FiveTuple fp = denspdf[id];
    FiveTuple fpp = denspdf[id+1];

    //Determine the probability distribution for galaxy density in this magnitude bin
    den_ent tpdf;
    for(int i=0;i<nbin;i++){
      tpdf.r[i] = (i+1)*rdel;
      float p1 = 0.5*(1. - fp[4])*(1+erf((log(tpdf.r[i])-fp[0])/(fp[1]*sqrt(2.0))));
      float p2 = 0.5*fp[4]*(1+erf((tpdf.r[i]-fp[2])/(fp[3]*sqrt(2.0))));
      float p3 = 0.5*(1. - fpp[4])*(1+erf((log(tpdf.r[i])-fpp[0])/(fpp[1]*sqrt(2.0))));
      float p4 = 0.5*fpp[4]*(1+erf((tpdf.r[i]-fpp[2])/(fpp[3]*sqrt(2.0))));
#ifdef SIXPARAMS
      p2 -= fpp[4]*2.0*t((tpdf.r[i]-fp[2])/fp[3],fp[5]);
      p4 -= fpp[4]*2.0*t((tpdf.r[i]-fpp[2])/fpp[3],fpp[5]);
#endif
      tpdf.prob[i] = weight1*(p1+p2) - weight2*(p3+p4);
    }

    //renormalize our probability distribution so that it integrates to unity
    float norm = tpdf.prob[nbin-1];
    if (norm > 0){
      for(int i=0;i<nbin;i++)
	tpdf.prob[i] /= norm;
    }
    pdf.push_back(tpdf);
  }
#ifdef DEBUG
  ofstream pdf_file("pdf_test.ascii");
  for(int id=0;id<nbin;id++){
    pdf_file<<x_prob[36].x[id]<<" "<<y_prob[36].x[id]<<endl;
  }
#endif

  cout<<endl;
}

/*
  Subroutine integrates our LF and gives us a list of galaxies
  with magnitude, local density, and redshift that needs to 
  be assigned to our particles
*/
vector <Galaxy *> GetGalaxies(double vol){

  //initialize the random number generator
  srand48(seed);

  float quasar_magmin = -18.0;
  float quasar_frac = 0.02;

  //setup factor for rescaling the LF
  //float rho = sim.Np()*sim.ParticleMass()/volume;
  //float rhoback = 3.*100*100/(8*!PI*4.301e-9)*sim.OmegaM();
  float phi_rescale = rho/rhoback;
  cout<<"Rescaling factor for LF normalization: "<<phi_rescale<<endl;


  //Setup our Luminosity function to be integrated
  integrator LF_integrator(&LF, 0, 1.0E-05, 0.01);

  //calculate the total number of galaxies we expect to generate
  int n = (int) (LumNumberDensity(Magmin)*vol);
  //do we need to add more galaxies for quasars to go onto?  
  if (Magmin < quasar_magmin){
    cout<<"Adding some extra galaxies for QSOs..."<<endl;
    int n1 = LumNumberDensity(Magmin)*vol*quasar_frac;
    int n2 = LumNumberDensity(quasar_magmin)*vol*quasar_frac;
    cout<<"  Added "<<n2-n1<<" new galaxies."<<endl;
    n += (n2-n1);
  }
  PRNT("GetGalaxies:", n);

  //Reserve the space for our galaxy vector
  vector <Galaxy*> galaxies;
  galaxies.reserve(n);
  double volume = vol;
  cout<<volume<<" "<<sim.CubeVolume()<<endl;

  //count the number of galaxies we've integrated so far
  int ngal = 0;

#ifdef DEBUG_GET_GALS
  //some sanity-checking files
  ofstream dens_file("dens_check.ascii"); //The magnitude, density, and redshift of each galaxy
  ofstream denspdf_file("denspdf_check.ascii"); //Checking the integration of our denspdf calculation
  ofstream gal_prop("galaxies.ascii"); //Just the magnitude and density of each of our galaxies
#endif

  //setup our denspdf distribution
  vector <FiveTuple> denspdf; //the 5 parameters that define our PDF
  vector <float> dmagbins; //our PDF is definied as a function of magnitude
  ReadPDFs(denspdf, dmagbins); 
  assert(denspdf.size()== dmagbins.size()); //sanity check

  //We're going to need to differentiate between bins
  //  bercause the PDF is a cumulative PDF for all galaxies
  //  brighter than dmabgins[i] -- setting up storage
  float pdfbinsize = dmagbins[0]-dmagbins[1]; //assume that the magnitudes are equally spaced
  int ind = dmagbins.size()-3; //Ignore the last few bins because we subtract between magnitudes (probably just needs to be one)
  FiveTuple fTup, fTup_prime;
  int weight1, weight2;
  float cur_mag = dmagbins[ind];
#ifdef DEBUG_GET_GALS
  cout<<"Checking the first denspdf values: ";
  denspdf[0].Print();
  cout<<" The pdf parameters (has "<<dmagbins.size()<<" elements)"<<endl;
#endif

  cout<<endl;

  //DENSPDF_IN_MAG_BIN should pretty mcuh always be set
#ifdef DENSPDF_IN_MAG_BIN
  //structures for storing our probability arrays: transform from cumulative to differential PDFs
  vector <den_ent> pdf;
  vector <int> weight;
  cout<<"Defining the probability arrays..."<<endl;
  define_prob(pdf, weight, denspdf, dmagbins, vol);
#ifdef DEBUG_GET_GALS
  for(int i=0;i<dmagbins.size();i++)
    denspdf_file<<weight[i]<<" "<<dmagbins[i]<<endl;
#endif
#endif

#ifdef SCALE_BY_GROWTH_FCN
  //finally, we need to know the growth
  vector <float> zarr;
  vector <float> darr;
  read_growthfcn2(zarr, darr);
#endif

  //no we can calculate magnitudes and densities
  cout<<"Generating Galaxy Magnitudes."<<endl;
  double ng_expected = LF_integrator.Integrate(-30, Magmin)*volume;
  cout<<"Expect to generate "<<ng_expected<<" galaxies."<<endl;

  //just set some initial values
  //  Note that ind was set earlier as dmagbins.size()-1:  
  //    the brightest objects we're interested in
  fTup = denspdf[ind];
  fTup_prime = denspdf[ind+1];
  weight1 = weight[ind];
  weight2 = weight[ind+1];
  cur_mag = dmagbins[ind];

  //generate galaxies by looping over our magnitude bins
  for(double mag=-23.7; mag <=-8; mag+=0.001){
    //test if it's time to advance to the next bin in our denspdf function
    if(mag>cur_mag){
      if(ind >= 0){
	if (ind != 0)  //note: we assume that the denspdf function is constant at dim magnitudes
	  --ind;
	fTup = denspdf[ind];
	fTup_prime = denspdf[ind+1];
#ifdef DENSPDF_IN_MAG_BIN
	weight1 = weight[ind];
	weight2 = weight[ind+1];
#endif
	cur_mag = dmagbins[ind];
	//cout<<ind<<" "<<mag<<" "<<cur_mag<<" "<<cur_mag+0.5*pdfbinsize<<" "<<galaxies.size()<<" "<<n<<endl;
      } else {
        fTup = denspdf[0];
        fTup_prime = denspdf[1];
        weight1 = 1.0;
        weight2 = 0.0;

	cout<<"[GetGalaxies] PDF binning problem... exiting"<<endl;
	cout<<ind<<" "<<mag<<" "<<cur_mag<<" "<<cur_mag+0.5*pdfbinsize<<endl;
	exit(1);
      }
    }

    //Determine the number of galaxies in this magnitude bin
    double num = LF_integrator.Integrate(-30, mag)*volume;
    //are we in the quasar regime?
    if (mag > Magmin) num *= quasar_frac;
    int diff = ((int) floor(num))-galaxies.size();

    //we allocated space for n galaxies -- truncating to stop if we're gone further than that
    if (num>n) diff = n-galaxies.size();

    //loop over new galaxies to determine redshifts and densities
    for(int i=0;i<diff;i++){
      float zGal = SelectGalaxyZ();

#ifdef DENSPDF_IN_MAG_BIN
      float dist8 = LocalDens(pdf[ind]);
      //have considered some models where we scale d8 by the growth function
#ifdef SCALE_BY_GROWTH_FCN
      float this_growth_fcn = 1.0;
      for (int i=0;i<zarr.size();i++){
        this_growth_fcn = darr[i];
        if (zarr[i] >= zGal)
          continue;
      }
      dist8 /= pow(this_growth_fcn,1.0/3.0);
#endif
      //we've defined our new galaxy:  save it and create it
      Galaxy * galaxy = new Galaxy(mag,ngal,dist8);
#else
      Galaxy * galaxy = new Galaxy(mag,ngal,fTup.LocalDens());
#endif
      galaxies.push_back(galaxy);
      galaxies[ngal]->zGal(zGal);
      ngal++;
    }
    if (galaxies.size()>=n) break;
  }  

  //output file with galaxy magnitudes, redshfits, and rdel's to confirm that code is working properly
#ifdef DEBUG_GET_GALS
  string galcheck_str = "galaxies_check.txt";
  ofstream galcheck_file(galcheck_str.c_str());
  for(int i=0;i<galaxies.size();i++)
    galcheck_file<<galaxies[i]->Mr()<<" "<<galaxies[i]->zGal()<<" "<<galaxies[i]->Dist8()<<endl;
#endif

  //output some diagnostics of our galaxy distribution
  cout<<"[GetGalaxies] Done getting galaxy magnitudes "<<galaxies.size()<<endl;
  int dimmest = -100.;
  int brightest = 0.;
  for(int ig=0;ig<galaxies.size();ig++){
    if (galaxies[ig]->Mr() < brightest)
      brightest = galaxies[ig]->Mr();
    if (galaxies[ig]->Mr() > dimmest)
      dimmest = galaxies[ig]->Mr();
  }
  cout<<"  brightest/dimmest galaxies = "<<brightest<<"/"<<dimmest<<endl;
  return galaxies;
}



double area(double alpha, double* dummy){
  return cos(alpha*3.141592654/180.);
}

//returns fractional area of the sky
double fractional_area(){
  
  integrator area_integrator(&area, 0,1.0e-03,0.01);
  double sky_coverage = (DECMAX-DECMIN)*area_integrator.Integrate(RAMIN,RAMAX);
  //  cout<<sky_coverage<<" "
  //  <<90.*area_integrator.Integrate(0., 90.)<<endl;
  return sky_coverage/41253.0;
}


