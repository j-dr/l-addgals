#include <vector>
#include <cassert> 
#include <unistd.h>     //for sleep()
#include "chainel.h"
#include "cosmo.h" 
#include "constants.h"
#include "myrand.h"  //for rand stuff
#include "stl_util.h"  //for copy if
#include "choose.h" //galaxy recipes
#include "galaxy.h"
#include "hv.h"


#include "sample.h"
#include "fivetuple.h"

using namespace std;
int xilike(string cffilename, string datfilename, double &chi2, int np);
int cf(string outcfn, string cffile);
int cf_fit(string outfilename, ofstream &datfile, double &r, double &g, double &chi);
int calculate_wp(string fname, string outfname, int nbins);

void DeleteAndNullifyBadgal(Galaxy*& pgalaxy);
//void AssignGalaxies(vector <Particle *>& particles, vector <Galaxy *>& galaxies);
void Assignment(vector <Particle *>& particles, vector <Galaxy *>& galaxies);
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

const  int nnbins = 50;
const  int mdbins = 24;


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

string mdlabel(string label, int iter){
  string runlab = MakeString(runcount*10+iter, 6);
  string dir = "nndat/"+label+"/";
  string s = dir+"meandens"+runlab+".dat";
  return s;
}

void calcnn(vector <Galaxy*> galaxies, string nnmfile, string mdmfile){
  vector <Galaxy *> galaxyslice;
  //note that you want to fill full zbins (=0.02)
  for(int gi=0;gi<galaxies.size();gi++){
    //    if((galaxies[gi]->Z()>0.06)&&(galaxies[gi]->Z()<=0.1))
    //    if((galaxies[gi]->Mr()<Magmin_col)&&(galaxies[gi]->Z()<=0.18))
    if((galaxies[gi]->Z()<=0.18))
    //    if((galaxies[gi]->Z()>0.00)&&(galaxies[gi]->Z()<=0.16))
      galaxyslice.push_back(galaxies[gi]);
    
  }
  vector <float> nndist = GetNeighborDist(galaxyslice);  
  //  cout<<galaxyslice.size()<<" "<<galaxies.size()
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

    double magmin = 22.2;
    vector <double> dens_mean(mdbins);
    vector <int> ninbin(mdbins);
    double chi2md;
    string mdfile = "meandens.dat";
    //string mdmfile = "outdat_"+label+"/meandens."+runlabel+".dat";
    ofstream mdmout(mdmfile.c_str());
    for(int di=0; di<mdbins; di++){
      dens_mean[di] = 0;
      ninbin[di] = 0;
    }
    for(int i = 0; i < galaxyslice.size(); i++){
      int magbin = (int) (floor((galaxies[i]->Mr()+magmin)/0.1));
      if(magbin<0) magbin=0;
      if((magbin>=0)&&(magbin<mdbins+2)){ 
	dens_mean[magbin] = dens_mean[magbin]+nndist[i];
	ninbin[magbin] = ninbin[magbin]+1;
      }
    }

    for(int di=0; di<mdbins; di++){
      //cout<<magmin*-1+di*0.1+0.05<<" "<<ninbin[di]<<endl;
      mdmout<<magmin*-1+di*0.1+0.05<<" "
	    <<dens_mean[di]/ninbin[di]<<" "
	    <<dens_mean[di]/ninbin[di]/sqrt(ninbin[di])<<endl;
    }
}

void calccf_external(vector <Galaxy*> galaxies, string label, int iter, Sample &sample){
  float *x, *y, *z;
  int thisn;
  string twp_file = './temp_gal_wpinfo.dat'
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
      twp_file<<p->X()<<" "<<p->Y()<<" "<<p->Z()<<" "<<p->Vx()<<" "<<p->Vy()<<" "<<p->Vz()<<endl;
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
      twp_file<<p->X()<<" "<<p->Y()<<" "<<p->Z()<<" "<<p->Vx()<<" "<<p->Vy()<<" "<<p->Vz()<<endl;
    }    
  }

  string dir = "cfdat/"+label+"/";
  string cffile = sample.Simfile(dir, runcount, iter);

  string cmd = "~/jeremy_wp/wp_covar "+sample.Rmin()+" "+sample.Rmax()+" "+sample.Nbins()+" "+sim.Boxszie()+" 0 "+sim.Boxsize()+" 1 "+twp_file+" a 1 1 1  > "+cffile;  
  system(cmd.c_str();
  string cmd = "rm "+twp_file;
  system(cmd.c_str();

  //covar(x, y, z, x, y, z, thisn, sim.Boxsize(), sample.Rmin(), sample.Rmax(), sample.Nbins(), cffile.c_str());
  if(sample.Cftype()==WP){
    string com = "more "+cffile;
    system(com.c_str());
    calculate_wp(cffile, cffile, sample.Nbins()-1);
    system(com.c_str());
  }
  //    else cout<<"sample"<<sample.Cftype()<<endl;
  free(x); free(y); free(z);
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
  extern int likemodel;
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
  Assignment(particles, galaxies);
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
    //    cout <<"chi2 "<<si<<": "<<chi2<<endl;
    params[si] = chi2;
  }
  if(likemodel == -1){
    chi2 = params[0]+params[1]+params[2];
    //    cout<<"using all cfs"<<params[0]<<" "<<params[1]<<" "<<params[2]<<" "<<chi2<<endl;
  }
  else if(likemodel == 0){
    chi2 = params[0];
  }
  else if(likemodel == 1)
    chi2 = params[1];
  else if(likemodel == 2)
    chi2 = params[2];
  else{
    string modellabel = nnlabel(label,iter);
    string modellabel2 = mdlabel(label,iter);
    //   string nnmfile = "outdat_"+label+"/nn5."+runlabel+".dat";
    calcnn(galaxies, modellabel, modellabel2);
    string nnfile = "nn5.dat";
    string mdfile = "meandens.dat";
    double chi2nn, chi2md;
    xilike(modellabel.c_str(), nnfile.c_str(), chi2nn, nnbins);
    params[3] = chi2nn;
    
    xilike(modellabel2.c_str(), mdfile.c_str(), chi2md, mdbins);
    params[4] = chi2md;
    if(likemodel == 3)
      //   chi2 = params[0]+params[1]+params[2]+params[3]+params[4];
      chi2 = params[0]+params[1]+params[2]+params[4];
  }
  /*    for(int si=0;si<samples.size();si++){
	string cffile = samples[si].Simfile(dir, runcount*10+iter);
	vector <float> xi_model_ave(samples[si].Nbins());
	  //      for(int iter=0;iter<niter;iter++){
	  
	  //    }
	  string cffile = samples[si].Simfile(dir, runcount*10+iter);
	  }

	  if(usenn){
	  string modellabel = nnlabel(label, iter);
	  xilike(nnmfile.c_str(), nnfile.c_str(), chi2nn, nnbins);
	  params[3] = chi2nn;
	  }
	  return chi2nn;
	  if(usenn){
	  params[3] = nnlike(galaxies, label, iter);

	  }
	  else 
	  
    */
}
