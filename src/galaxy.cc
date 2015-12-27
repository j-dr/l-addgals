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
#include "galaxy_global.h"
#ifdef LF_FROM_DATA
#include <vector>
#endif

using namespace std;

//NBIN is the total number of rdel bins for our magnitude calculation
//#define NBIN 85000
//#define NBIN 8500
//struct den_ent{
//  float prob[NBIN];
//  float r[NBIN];
//};

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
  float rho = sim.Np()*sim.ParticleMass()/vol;
  float rhoback = 3.*100*100/(8*!PI*4.301e-9)*sim.OmegaM();
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
  //structures for storing our probability arrays: transform from cumulative to differential PDFs
  vector <den_ent> pdf;
  vector <int> weight;
  cout<<"Defining the probability arrays..."<<endl;
  define_prob(pdf, weight, denspdf, dmagbins, vol);
#ifdef DEBUG_GET_GALS
  for(int i=0;i<dmagbins.size();i++)
    denspdf_file<<weight[i]<<" "<<dmagbins[i]<<endl;
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

