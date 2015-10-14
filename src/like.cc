/*#include <vector>
#include <cassert> 
#include <unistd.h>     //for sleep()

#include "cosmo.h" 
#include "constants.h"
#include "myrand.h"  //for rand stuff
#include "stl_util.h"  //for copy if
#include "choose.h" //galaxy recipes
//#include "galaxy.h"
#include "hv.h"

#include "chainel.h"
*/


using namespace std;

int main(){
  int a = 1;
  -11-:---F1  like.cc        Top L15     (C++ Abbrev)-----------------------------

File Edit Options Buffers Tools C++ Help                                        
#include "cosmo.h"                                                              
#include "constants.h"                                                          
#include "myrand.h"  //for rand stuff                                           
#include "stl_util.h"  //for copy if                                            
#include "choose.h" //galaxy recipes                                            
    //#include "galaxy.h"                                                           
#include "hv.h"                                                                 
                                                                                
#include "chainel.h"                                                            
*/
#include "sample.h"

#include "fivetuple.h"

    using namespace std;

  int main(){
    int a = 1;
#include "cosmo.h"                                                              
#include "constants.h"                                                          
#include "myrand.h"  //for rand stuff                                           
#include "stl_util.h"  //for copy if                                            
#include "choose.h" //galaxy recipes                                            
      //#include "galaxy.h"                                                           

                                                                                

*/
      //#include "sample.h"
      //#include "chainel.h"                                                            

#include "hv.h"                                                                 #include "fivetuple.h"
using namespace std;

int main(){
  int a = 1;
  a = 1+1;
  //  cout<<"hello"<<endl;

}


/*int xilike(string cffilename, string datfilename, double &chi2, int np);
int cf(string outcfn, string cffile);
int cf_fit(string outfilename, ofstream &datfile, double &r, double &g, double &chi);
int calculate_wp(string fname, string outfname, int nbins);

void DeleteAndNullifyBadgal(Galaxy*& pgalaxy);
void AssignGalaxies(vector <Particle *>& particles, vector <Galaxy *>& galaxies);
void UnAssignGalaxies(vector <Particle *>& particles, vector <Galaxy *>& galaxies);
void cfcovar(string infn, string outfn){
  string command = "covar3 0.127 20 12 375.0 0 375.0 1 "+infn+" a 1 1 auto | awk '{if (NR>2) print $3, $4, $5}' > "+outfn;
  system(command.c_str());
}

vector <float> GetNeighborDist(vector <Galaxy *> galaxies);   

extern "C" void covar(float *x1, float *y1, float *z1, 
		      float *vx1, float *vy1, float *vz1, int np1,
		      float rcube, float rmin, float rmax, int nbin, const char *filename);


extern "C" void correlation(float *x1, float *y1, float *z1, int np1, 
		      float rcube, float rmin, float rmax, int nbin, const char *filename);

void Delete(Galaxy*& pgalaxy){
    delete pgalaxy;
    pgalaxy =0;
}

static int nbri = 0;
static int nmid = 0;
static int ndim = 0;
static int runcount = 0;

static int midpoint;
static int dimpoint;

// compare galaxy pointers by local mass density
bool BrighterGal(Galaxy * a, Galaxy * b)
{
  return a->Mr() < b->Mr(); 
}

string nnlabel(string label, int iter){
  string runlab = MakeString(runcount*10+iter, 6);
  string dir = "nndat/"+label+"/";
  string s = dir+"nn5"+runlab+".dat";
  return s;
}

void calcnn(vector <Galaxy*> galaxies, string nnmfile){
  vector <Galaxy *> galaxyslice;
  //note that you want to fill full zbins (=0.02)
  for(int gi=0;gi<galaxies.size();gi++){
    if((galaxies[gi]->Z()>0.06)&&(galaxies[gi]->Z()<=0.1))
      galaxyslice.push_back(galaxies[gi]);
  }
  vector <float> nndist = GetNeighborDist(galaxyslice);  
  //  cout<<galaxyslice.size()<<" "<<galaxies.size()
  string nnfile = "nn5.dat";

  int nnbins = 50;
  float binsize=0.1;
    vector <float> nnval(nnbins);
    vector <float> nnerr(nnbins);
    for(int i=0; i<nnbins;i++){
      nnval[i] = 0;
    }
    for(int gi=0; gi<nndist.size();gi++){
      int this_bin = (int) (floor((nndist[gi]-0.5*binsize)/binsize));
      if ((this_bin>=0)&&(this_bin<nnbins)){
	nnval[this_bin] = nnval[this_bin]+1;
      }
    }
    ofstream nnout(nnmfile.c_str());
    for(int i=0; i<nnbins;i++){
      float tmp = nnval[i];
      nnval[i] = tmp/nndist.size();
      nnerr[i] = sqrt(tmp)/nndist.size();
      nnout<<(i+1)*binsize<<" "<<nnval[i]<<" "<<nnerr[i]<<endl;
    }
}

void calccf(vector <Galaxy*> galaxies, string label, int iter, Sample &sample){
  float *x, *y, *z;
  int thisn;
  if(sample.Ncf()<sample.Nbin()){
    thisn = sample.Ncf();
    x=(float *)malloc(thisn*sizeof(float)) ;
    y=(float *)malloc(thisn*sizeof(float)) ;
    z=(float *)malloc(thisn*sizeof(float)) ;
    for(int gi=0;gi<thisn;gi++){
      int midi = randint(sample.Startpoint(),sample.Endpoint());
      Galaxy* gal = galaxies[midi];
      Particle* p = gal->P();
      x[gi] = p->X();
      y[gi] = p->Y();
      z[gi] = p->Z();
    }          
  }
    else{
    thisn = sample.Ncf();
      x=(float *)malloc(thisn*sizeof(float)) ;
      y=(float *)malloc(thisn*sizeof(float)) ;
      z=(float *)malloc(thisn*sizeof(float)) ;
      for(int gi=0;gi<thisn;gi++){
	Galaxy* gal = galaxies[gi];
	Particle* p = gal->P();
	x[gi] = p->X();
	y[gi] = p->Y();
	z[gi] = p->Z();
      }    
    }
    string dir = "cfdat/"+label+"/";
    string cffile = sample.Simfile(dir, runcount, iter);
    covar(x, y, z, x, y, z, thisn, sim.Boxsize(), sample.Rmin(), sample.Rmax(), sample.Nbins(), cffile.c_str());
    if(sample.Cftype()==WP){
      string com = "more "+cffile;
      system(com.c_str());
      calculate_wp(cffile, cffile, sample.Nbins()-1);
      system(com.c_str());
    }
    //    else cout<<"sample"<<sample.Cftype()<<endl;
    free(x); free(y); free(z);
 }

void CountGals(vector <Galaxy *> galaxies){
  if(ndim==0){
    for(int gi=0;gi<galaxies.size();gi++){
      if((galaxies[gi]->Mr()>=-23)&&(galaxies[gi]->Mr()<-22))
	nbri++;
      else if((galaxies[gi]->Mr()>=-22)&&(galaxies[gi]->Mr()<-21)){
	if(nmid == 0) midpoint = gi;
	nmid++;
      }
      else if((galaxies[gi]->Mr()>=-21)&&(galaxies[gi]->Mr()<-20)){
	if(ndim == 0) dimpoint=gi;
	ndim++;
      }
    }
  }
}

void Likelihood(vector <Particle *> particles, vector <Galaxy *> galaxies, 
		ChainEl chel,
		vector <double> &params,
		double &chi2, string label){
  runcount++;
  //  int niter = 2;
  bool usenn = false;

  CountGals(galaxies);
  vector <Sample> samples;
  sample_type samplet = BRIGHT;
  cf_type cft = XI;
  Sample sampleb(samplet, nbri, 0, midpoint, cft);

  samplet = MID;
  Sample samplem(samplet, nmid, midpoint, dimpoint,cft);
  
  samplet = DIM;
  Sample sampled(samplet, ndim, dimpoint, galaxies.size()-1,cft);

  samples.push_back(sampleb);
  samples.push_back(samplem);
  samples.push_back(sampled);
  int iter = 0;
  //  for(int iter=0;iter<niter;iter++){

  for(int gi=0;gi<galaxies.size();gi++){
    float mr = pow(10.,-0.4*(galaxies[gi]->Mr()+20.44));
    FiveTuple fTup(chel.cmean(),chel.fmean(mr),chel.fsig(mr),chel.ffrac(mr)); 
    galaxies[gi]->Dist8(fTup.LocalDens());
  }
  AssignGalaxies(particles, galaxies);
  // assign sorts; need to sort this back
  sort(galaxies.begin(),galaxies.end(),BrighterGal);
  //  cout<<"ss"<<samples.size()<<endl;
  for(int si=0;si<samples.size();si++){
    int snb = samples[si].Nbins();
    //    cout<<si<<" "<<snb<<" "<<label<<" "<<iter<<endl;
    calccf(galaxies, label, iter, samples[si]);
    string dir = "cfdat/"+label+"/";
    string cffile = samples[si].Simfile(dir, runcount, iter);
    ifstream cfinfile(cffile.c_str());
    vector <float> rr(snb);
    vector <float> xi(snb);
    vector <float> xierr(snb);

    for(int bi=0;bi<snb-2;bi++){
      cfinfile>>rr[bi]>>xi[bi]>>xierr[bi];
      //      cout<<"h"<<rr[bi]<<" "<<xi[bi]<<endl;
    }
    cffile = samples[si].Simfile(dir, runcount);
    ofstream outfile(cffile.c_str());
    for(int bi=0;bi<samples[si].Nbins()-2;bi++){
      outfile<<rr[bi]<<" "<<xi[bi]<<" "<<xierr[bi]<<endl;
    }
    string datfile = samples[si].Datfile();
    xilike(cffile.c_str(), datfile.c_str(), chi2, samples[si].Nbin());
    //cout <<"chi2"<<chi2<<endl;
    params[si] = chi2;
  }


  if(usenn){
   string modellabel = nnlabel(label,iter);
   calcnn(galaxies, modellabel);
  }
  
  //}

  //chi2 = params[0]+params[1]+params[2];
  chi2 = params[0];
}
*/
