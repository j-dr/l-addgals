#include <cmath>
#include <iostream>
#include "nr.h"
#include "myrand.h"
#include "fivetuple.h"
#include "constants.h"
#include "galaxy.h"
#include "cosmo.h"
#include "ReadParameters.h"
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


//get the explicit denpdf parameters for our redshift, magnitude
FiveTuple denspdf_params(float magnitude, float zred)
{

  float mag_ref= -20.5;
  //float dim_mag_lim = -19.0-mag_ref;
  float bright_mag_lim = -22.5-mag_ref;
  float dim_mag_lim  = -18 - mag_ref;
  float mag = magnitude-mag_ref;
  
  //redshift evolution fit in terms of scale factor shifted 
  //to make pivot close to mean scale factor of BCC boxes

  float a   = 1 / (zred + 1) - 0.35;
  if (mag > dim_mag_lim) mag = dim_mag_lim;
  if (mag < bright_mag_lim) mag = bright_mag_lim;
  float cm, cs, fm, fs, p;

  cm = cm0 + cm1*mag + cm2*mag*mag + cm3*mag*mag*mag + cm4*mag*mag*mag*mag + cmz1*a + cmz2*a*a
         + cmz3*a*a*a + cmz4*a*a*a*a + cm1z1*mag*a + cm1z2*mag*a*a + cm2z1*mag*mag*a 
         + cm2z2*mag*mag*a*a + cm3z1*mag*mag*mag*a + cm1z3*mag*a*a*a;
  cs = cs0 + cs1*mag + cs2*mag*mag + cs3*mag*mag*mag + cs4*mag*mag*mag*mag + csz1*a + csz2*a*a
         + csz3*a*a*a + csz4*a*a*a*a + cs1z1*mag*a + cs1z2*mag*a*a + cs2z1*mag*mag*a 
         + cs2z2*mag*mag*a*a + cs3z1*mag*mag*mag*a + cs1z3*mag*a*a*a;
  fm = fm0 + fm1*mag + fm2*mag*mag + fm3*mag*mag*mag + fm4*mag*mag*mag*mag + fmz1*a + fmz2*a*a
         + fmz3*a*a*a + fmz4*a*a*a*a + fm1z1*mag*a + fm1z2*mag*a*a + fm2z1*mag*mag*a 
         + fm2z2*mag*mag*a*a + fm3z1*mag*mag*mag*a + fm1z3*mag*a*a*a;
  fs = fs0 + fs1*mag + fs2*mag*mag + fs3*mag*mag*mag + fs4*mag*mag*mag*mag + fsz1*a + fsz2*a*a
         + fsz3*a*a*a + fsz4*a*a*a*a + fs1z1*mag*a + fs1z2*mag*a*a + fs2z1*mag*mag*a 
         + fs2z2*mag*mag*a*a + fs3z1*mag*mag*mag*a + fs1z3*mag*a*a*a;
  p = p0 + p1*mag + p2*mag*mag + p3*mag*mag*mag + p4*mag*mag*mag*mag + pz1*a + pz2*a*a
         + pz3*a*a*a + pz4*a*a*a*a + p1z1*mag*a + p1z2*mag*a*a + p2z1*mag*mag*a 
         + p2z2*mag*mag*a*a + p3z1*mag*mag*mag*a + p1z3*mag*a*a*a;

  if (p<0) p=0;
  if (p>1) p=1;
  if (cs<0) cs=0;
  if (fm<0) fm=0;
  if (fs<0) fs=0;

  FiveTuple fp(cm,cs,fm,fs,p);

  return fp;
}



/*
  Subroutine that returns the Probability for a galaxy with magnitude M to have local density rdel
  by converting from our magnitude limited denspdf to one in a magnitude bin.  Requires information
  on the global luminosity function
*/
den_ent define_prob(float magnitude, float redshift, float vol)
{

  //cout<<"Starting define_prob..."<<endl;
  float dmag = 0.05;

  //Note that NBIN defines the total number of probability bins we'll be using in r
  int nbin = NBIN;
  float rdel = 8.5/nbin;

  //The difference in weight1 and weight2 gives the # of galaxies in this mag bin
  //cout<<"Getting weights..."<<endl;
  float weight1 = (int) (LumNumberDensity(magnitude+dmag)*vol);
  float weight2 = (int) (LumNumberDensity(magnitude-dmag)*vol);
  //cout<<"weight1 = "<<weight1<<", weight2 = "<<weight2<<endl;

  //The CUMULATIVE density PDFs both above and below our current magnitude bin
  //cout<<"Getting the denspdf parameters..."<<endl;
  FiveTuple fp = denspdf_params(magnitude+dmag, redshift);
  FiveTuple fpp = denspdf_params(magnitude-dmag, redshift);
#ifdef DEBUG_DENSPDF
  ofstream pdf_file("denspdf_param_test.ascii", std::ofstream::out | std::ofstream::app);
  pdf_file<<magnitude+dmag << " "<< redshift << " "<< fp[0] << " "<< fp[1] << " "<< fp[2] << " "<< fp[3] << " "<< fp[4] << " "<< endl;
#endif
  //if (magnitude+dmag > -19.0)
  //  fp = denspdf_params(-19.0, redshift);
  //if (magnitude-dmag > -19.0)
  //  fpp = denspdf_params(-19.0, redshift);
  //if (magnitude > -19.0) {
  //  fp = denspdf_params(-19.0+dmag, redshift);
  //  fpp = denspdf_params(-19.0-dmag, redshift);
  //}
  //cout<<"Parameters 1: "<<fp[0]<<", "<<fp[1]<<", "<<fp[2]<<", "<<fp[3]<<", "<<fp[4]<<endl;
  //cout<<"Parameters 2: "<<fpp[0]<<", "<<fpp[1]<<", "<<fpp[2]<<", "<<fpp[3]<<", "<<fpp[4]<<endl;

  //Determine the probability distribution for galaxy density in this magnitude bin
  //cout<<"Making the probability distribution..."<<endl;
  den_ent pdf;
  for(int i=0;i<nbin;i++){
    pdf.r[i] = (i+1)*rdel;
    float p1 = 0.5*(1. - fp[4])*(1+erf((log(pdf.r[i])-fp[0])/(fp[1]*sqrt(2.0))));
    float p2 = 0.5*fp[4]*(1+erf((pdf.r[i]-fp[2])/(fp[3]*sqrt(2.0))));
    float p3 = 0.5*(1. - fpp[4])*(1+erf((log(pdf.r[i])-fpp[0])/(fpp[1]*sqrt(2.0))));
    float p4 = 0.5*fpp[4]*(1+erf((pdf.r[i]-fpp[2])/(fpp[3]*sqrt(2.0))));
    pdf.prob[i] = weight1*(p1+p2) - weight2*(p3+p4);
  }

  //renormalize our probability distribution so that it integrates to unity
  //cout<<"Normalizing the distribution..."<<endl;
  float norm = pdf.prob[nbin-1];
  if (norm > 0){
    for(int i=0;i<nbin;i++)
      pdf.prob[i] /= norm;
  }
#ifdef DEBUG
  ofstream pdf_file("pdf_test.ascii");
  for(int id=0;id<nbin;id++){
    pdf_file<<x_prob[36].x[id]<<" "<<y_prob[36].x[id]<<endl;
  }
#endif

  return pdf;
}


float LocalDens(den_ent pdf)
{


  float d8_max = 8.5;
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
      ind--;
    float dx = ranu - pdf.prob[ind-1];
    float slope = (pdf.r[ind]-pdf.r[ind-1])/(pdf.prob[ind]-pdf.prob[ind-1]);
    d8 = pdf.r[ind-1]+dx*slope;
  }
  if(!(d8 < d8_max)) d8 = d8_max;
  if(!(d8 > d8_min)) d8 = d8_min;
  return d8;
}


float LocalDens(float magnitude, float redshift, float vol)
{

  //define our probability distribution for this galaxy
  //cout<<"Defining probability array..."<<endl;
  den_ent pdf = define_prob(magnitude, redshift, vol);
  //for(int i=0;i<NBIN;i++){
    //cout<<pdf.r[i]<<" "<<pdf.prob[i]<<endl;
  //}

  float d8 = LocalDens(pdf);

  return d8;
}


float SelectGalaxyZ()
{
  //selects a random radius for galaxy s.t. prob \propto r^3
#ifdef SNAPSHOT
  float z = sim_redshift;
#else
  float rn = drand48();
  float rm = cosmo.RofZ(ZREDMAX);
  float rmin = cosmo.RofZ(ZREDMIN);
  //rn *= rm*rm*rm;
  rn = (rm*rm*rm - rmin*rmin*rmin)*rn + rmin*rmin*rmin;
  rn = pow(rn, 1.0/3.0);
  float z = cosmo.ZofR(rn);
#endif
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

vector <float> gmagbin;
vector <float> gdensity;
#ifdef LF_FROM_DATA

void read_lf_data(void){
  string lf_file = "LF.dat";
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
    gdensity.push_back(gdensity1);
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
  Subroutine integrates our LF and gives us a list of galaxies
  with magnitude, local density, and redshift that needs to
  be assigned to our particles
*/
vector <Galaxy *> GetGalaxies(double vol, float phi_rescale){

  //initialize the random number generator
  srand48(seed);


  //rescale the LF for cosmic variance
  for(int i=0;i<gdensity.size();i++)
    gdensity[i] *= phi_rescale;

  //Setup our Luminosity function to be integrated
  integrator LF_integrator(&LF, 0, 1.0E-05, 0.01);

  //calculate the total number of galaxies we expect to generate
  int n = (int) (LumNumberDensity(Magmin)*vol);
  PRNT("GetGalaxies:", n);
  //cout<<"Simulation redshift: "<<sim_redshift<<endl;

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

#ifdef SCALE_BY_GROWTH_FCN
  //finally, we need to know the growth
  vector <float> zarr;
  vector <float> darr;
  read_growthfcn2(zarr, darr);
#endif

  //create an array of pdf's for rapid lookup
  float dz_pdf = 0.01;
  int nzbins = 0;
  for (float tz = ZREDMIN; tz <= ZREDMAX; tz+=dz_pdf)
    nzbins++;

  //no we can calculate magnitudes and densities
  cout<<"Generating Galaxy Magnitudes."<<endl;
  double ng_expected = LF_integrator.Integrate(-30, Magmin)*volume;
  cout<<"Expect to generate "<<ng_expected<<" galaxies."<<endl;
  float next = -23.0;

  //generate galaxies by looping over our magnitude bins
  for(double mag=-25; mag <=-8; mag+=0.001){
    if (mag > next){
      cout<<"Generating galaxy magnitudes "<<mag<<endl;
      next += 1;
    }

    //fill an array of pdf's for rapid lookup
    vector <den_ent> pdfarray;
    float tz = 0.0;
    float tmag;
    pdfarray.reserve(nzbins);
    for (int iz=0;iz<nzbins;iz++)
      {
	tz = ZREDMIN + iz*dz_pdf;
	tmag = mag;
	if (tmag > -19.0) tmag = -19.0;
        pdfarray[iz] = define_prob(tmag,tz,vol);
      }


    //Determine the number of galaxies in this magnitude bin
    double num = LF_integrator.Integrate(-30, mag)*volume;
    int diff = ((int) floor(num))-galaxies.size();
    //cout<<" Will make "<<diff<<" galaxies with magnitude "<<mag<<endl;

    //we allocated space for n galaxies -- truncating to stop if we're gone further than that
    if (num>n) diff = n-galaxies.size();

    //loop over new galaxies to determine redshifts and densities
    for(int i=0;i<diff;i++){
      //cout<<"Getting redshift..."<<endl;
      float zGal = SelectGalaxyZ();

      //cout<<"Getting local density with mag = "<<mag<<", zGal = "<<zGal<<", vol = "<<vol<<endl;
      //float dist8 = LocalDens(mag, zGal, vol);
      int pdf_zbin = round((zGal-ZREDMIN)/dz_pdf);
      if (pdf_zbin < 0) pdf_zbin = 0;
      if (pdf_zbin >= nzbins) pdf_zbin = nzbins-1;
      //cout<<"Getting local density with z = "<<zGal<<" and pdf_zbin = "<<pdf_zbin<<", corresponding to "<<ZREDMIN+dz_pdf*pdf_zbin<<endl;
      float dist8 = LocalDens(pdfarray[pdf_zbin]);

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
      //cout<<"Saving Galaxy..."<<endl;
      Galaxy * galaxy = new Galaxy(mag,ngal,dist8);
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


void read_out_galaxy_info(vector<Galaxy *> &gal, vector<float> &mr,
			  vector<float> &z, vector<GalSED> &seds,
			  vector<int> &sed_id, vector<int> &sed_cat_id)
{
  int i;
  vector<Galaxy *>::iterator itr;
  Particle * p;
  for (i=0; i<gal.size(); i++)
    {
      p = gal[i]->P();
      mr[i] = gal[i]->Mr();
      z[i] = p->ZredReal();
      sed_cat_id[i] = seds[sed_id[i]].CatId();
    }
}

void read_out_galaxy_info_w_densities(vector<Galaxy *> &gal, vector<float> &mr, vector<float> &z,
				      vector<GalSED> &seds, vector<int> &sed_id, vector<int> &sed_cat_id,
				      vector<float> &dist8)
{
  int i;
  vector<Galaxy *>::iterator itr;
  Particle * p;
  for (i=0; i<gal.size(); i++)
    {
      p = gal[i]->P();
      mr[i] = gal[i]->Mr();
      z[i] = p->ZredReal();
      sed_cat_id[i] = seds[sed_id[i]].CatId();
      dist8[i] = gal[i]->Dist8();
    }
}
