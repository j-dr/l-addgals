
#include <vector>
#include <cassert> 
#include <unistd.h>     //for sleep()
#include "nr.h"
#include "cosmo.h" 
#include "constants.h"
#include "myrand.h"  //for rand stuff
#include "stl_util.h"  //for copy if
#include "choose.h" //galaxy recipes
#include "hv.h"


/*
#include <vector>
#include <cassert> 
#include <unistd.h>     //for sleep()
#include <fstream>
#include "hv2.h"
#include "hv.h"
#include "fivetuple.h"
*/



//#define LOOP
//#define LIST

#include "chainel.h"
#include <iomanip.h>
#define MCMC
#define HOD
using namespace std;

#ifdef HVL
Simulation sim(box, use_cells);
#endif
#ifdef LANL
Simulation sim(box);
#endif
#ifdef GADGET
Simulation sim(box);
#endif
#ifdef MS_Gas
Simulation sim(box);
#endif
#ifdef CarmenLC
Simulation sim(box);
#endif


void CalculateCovar(vector <ChainEl> chain, string covfilename, int skipfirst);

/*
double r0ofmag(double mag){
  const double r0p[4] = {-1810.73, -282.367, -14.6294, -0.252483};
  return r0p[0]+r0p[1]*mag+r0p[2]*mag*mag+r0p[3]*mag*mag*mag;
}

double gamofmag(double mag){
  const double gamp[3] = {1.87398,-0.0386048,0.0135761};
  double lum =  pow(10.,(-0.4*(mag+20.44)));
  return -1*(gamp[0]+gamp[1]*lum+gamp[2]*lum*lum);
}
*/

Cosmology cosmo = sim.SimCosmology();

vector <Particle *> ReadParticles(int &nread);

bool DLess(Particle * a, Particle * b)
{
  return a->Dist8() < b->Dist8(); 
}

/*
void pdf_to_xi(vector <Particle *> particles, vector <Galaxy *> galaxies,
	       ChainEl chainel,
	       //double &r0, double &gam, 
	       vector <double> &params,
	       double &chi2, ofstream &datfile, string label);
*/
void Likelihood(vector <Particle *> particles, vector <Galaxy *> galaxies,
	       ChainEl chainel,
	       vector <double> &params,
	       double &chi2, string label);

static int runcount;
int likemodel;
int main(int argc, char**argv){
  const int Nchain =101000;
  cosmo.Print();
  cout<<sim.Boxsize()<<endl;
  cout<<sim.Boxsize()*BOXFR<<endl;
  cosmo.ReadZFile();
  if((argc<2)||(argc>3)){
    cout<<"Usage:pdfxi filelabel"<<endl;
    cout<<"or Usage:pdfxi filelabel likelihood model"<<endl;
    exit(1);
  }
  else if(argc==3){
    likemodel = atoi(argv[2]);
    cout<<"using likelihood model"<<likemodel<<endl;
  }
  else{
    likemodel = -1;
    cout<<"using likelihood model"<<likemodel<<endl;
  }
  string label = argv[1];
  string outfn = "./output_"+label+"_31k.dat";
  string com1 = "mkdir cfdat/"+label+"/";
  string com2 = "mkdir nndat/"+label+"/";
  system(com1.c_str());
  system(com2.c_str());
  // Read in particles from cubes.
  int nread =0;
  vector <Particle*> particles = ReadParticles(nread);
  float volume_fraction = particles.size()*1.0/nread;
  cout<<" Read "<<particles.size()<<" particles"<<endl;
  assert(particles.size()>0);

  cout<<"Sorting"<<endl;
  // Sort particles on local mass density
  sort(particles.begin(),particles.end(),DLess);
  vector <Galaxy *> galaxies;     
  //  galaxies=GetGalaxiesNodens(sim.CubeVolume()*volume_fraction);   
  //Calculate the volume of the region.  Convert to radians and switch dec to standard spherical units
  float PI = 3.1415926535897932384626433832795;
  float decmin_rad = PI*(90 - DECMAX)/180;
  float decmax_rad = PI*(90 - DECMIN)/180;
  float ramin_rad = PI*RAMIN/180;
  float ramax_rad = PI*RAMAX/180;
  float dra = ramax_rad - ramin_rad;
  float r_zmin = cosmo.RofZ(ZREDMIN);
  float r_zmax = cosmo.RofZ(ZREDMAX);
  float volume = (ramax_rad-ramin_rad)*(cos(decmin_rad)-cos(decmax_rad))*(pow(r_zmax,3)-pow(r_zmin,3))/3.;
#ifdef SNAPSHOT
  volume = sim.LengthUnit()*sim.LengthUnit()*sim.LengthUnit();
#endif
  cout<<"Volume in lightcone: "<<volume<<endl;
  //  galaxies=GetGalaxies(volume);
  galaxies=GetGalaxiesNodens(volume);   

  float acctot = 0.;
  double oldlik=0, newlik=0, accept=0;
  ChainEl oldel, newel;
  string chain_out = "chain_mcmc"+label+"_"+simlabel+".dat";
  string chi_out = "chi2_mcmc_"+label+"_"+simlabel+".dat";
  string covfile = "steps_"+label+"_"+simlabel+".dat";
  ofstream chout(chain_out.c_str());
  ofstream chiout(chi_out.c_str());

  vector <ChainEl> chain;
  int skipfirst = 0;
  for(int i=0; i<Nchain; i++){
    runcount++;
    bool acc = false;
    double chi2 = 10000;
    vector <double> meanparams(5);
    vector <double> params(5);
    chain.push_back(newel);
    if((chain.size()%1000==0)&&(chain.size()<Nchain/2)){
      if(chain.size()==2000)  skipfirst = 1000;
      //newel.PrintRotation();
      CalculateCovar(chain, covfile, skipfirst);
      //newel.PrintRotation();
    }
    if(newel.Okay()==-1){
      //int niter = 2;
      //for(int iter=0;iter<niter;iter++){
      Likelihood(particles, galaxies, newel, params, chi2, label);
      //cout<<iter<<" ";
      //for(int ii=0;ii<3;ii++){
      //cout<<params[ii]<<" ";
      //meanparams[ii]+=params[ii];
      //}
      //cout<<endl;
      //}
      //for(int ii=0;ii<3;ii++)
      //meanparams[ii]/=niter;
      //chi2 = meanparams[0]+meanparams[1]+meanparams[2];

      accept = 0;
      if(i==0){
	oldlik = -0.5*chi2;
      }
      else{
	newlik = -0.5*chi2;
	if (newlik > oldlik){
	  accept = 1.0;
	  acc = true;
	}
	else
	  accept = exp(newlik-oldlik);
      }
      if(NR::ran1(seed)<accept){
	acc = true;
      }
    }
    else{
      acc = 0;
    }
    //if(i==0) acc = 1;
    if(i==0)chout<<i<<" "<<1<<" "; else chout<<i<<" "<<acc<<" ";
    if(acc) newel.Write(chout);
    else oldel.Write(chout);
    chout<<endl;
    //    for(int jj=0;jj<4;jj++) chout<<meanparams[jj]<<" ";
    //write the chi_squared information
    for(int jj=0;jj<params.size();jj++) 
      chiout<<params[jj]<<" ";
    chiout<<chi2<<" "<<oldlik<<" "<<newlik<<" "<<newel.Okay()<<endl;

    acctot +=acc;
    cout<<i<<" "<<acc<<" ";
    newel.Print();
    cout<<acctot/(i+1.)<<" "<<newel.Okay()<<" "<<oldlik<<" "<<newlik<<" "<<chi2<<endl;
    if (acc) {
      oldel = newel;
      oldlik = newlik;
    }
    newel.Reset(oldel);
    for(int pi=0; pi<particles.size();pi++){
      particles[pi]->UnMakeGal();
    }
  }
#ifdef WRITECF
  string com = "cp "+cffile+" "+"cfdat_"+label+"/cf"+MakeString(i,4)+".dat";
  system(com.c_str());
#endif
  //  i++;
  cout<<"done with program."<<endl;
}

