#include <vector>
#include <cassert> 
#include "cosmo.h" 
#include "constants.h"
#include "myrand.h"  //for rand stuff
#include "stl_util.h"  //for copy if
#include <unistd.h>     //for sleep()
#include "choose.h" //galaxy recipes
#include "hv.h"
//#include "/home/risa/code/recipes/nr.h"
#include "/afs/slac.stanford.edu/u/ki/mbusha/projects/addgals/RisaLibs/recipes/nr.h"
//#define LOOP
//#define LIST
#define MCMC
#define HOD
using namespace std;

#ifdef HVL
Simulation sim(box, use_cells);
#endif
#ifdef LANL
Simulation sim(box);
#endif

double r0ofmag(double mag){
  const double r0p[4] = {-1810.73, -282.367, -14.6294, -0.252483};
  return r0p[0]+r0p[1]*mag+r0p[2]*mag*mag+r0p[3]*mag*mag*mag;
}

double gamofmag(double mag){
  const double gamp[3] = {1.87398,-0.0386048,0.0135761};
  double lum =  pow(10.,(-0.4*(mag+20.44)));
  return -1*(gamp[0]+gamp[1]*lum+gamp[2]*lum*lum);
}


Cosmology cosmo = sim.SimCosmology();

vector <Particle*> ReadParticles(void);

bool DLess(Particle * a, Particle * b)
{
  return a->Dist8() < b->Dist8(); 
}

void pdf_to_xi(vector <Particle *> particles, 
	       double cmean, double csig, double fmean, double fsig, double ffrac, 
	       double &r0, double &gam, double &chi2, ofstream &datfile, string label);



#ifdef MCMC
static const float sig[5]={0.34, 0.11, 0.39, 0.50, 0.22};

#define	GRND	sqrt(-2.*log(drand48()))*cos(2*M_PI*drand48())
class ChainEl{

 public:
  //start with brightest galaxy parameters
  ChainEl(){
    cmean = -0.4;
    csig = 0.46;
    fmean = 0.9;
    fsig = 1.47;
    ffrac = 0.5;
    }

  float cmean;
  //note::"sig" is actually the variance.
  float csig;
  float fmean;
  float fsig;
  float ffrac;

  void Reset(ChainEl oldel){
    cmean = oldel.cmean;
    //    cmean = oldel.cmean+ sig[0]*GRND ;
    //if(cmean<-1) cmean = -1;
    //if(cmean>0) cmean = 0;
    csig  =  oldel.csig;
      //+ sig[1]*GRND  ;
      //if(csig<0.3) csig = 0.3;
      //if(csig>0.5) csig = 0.5;
    fmean =  oldel.fmean + sig[2]*GRND  ;
    if(fmean<0.1) fmean = 0.1;
    if(fmean>4.0) fmean = 4.0;
    //if(fmean>2.0) fmean = 2.0;
    fsig  =  oldel.fsig + sig[3]*GRND  ;
    if(fsig<0.2) fsig = 0.2;
    //if(fsig>2.2) fsig = 2.2;
    if(fsig>3.0) fsig = 3.0;
    ffrac =  oldel.ffrac + sig[4]*GRND  ;
    if(ffrac<0.05) ffrac = 0.05;
    if(ffrac>1.0) ffrac = 1.0;
  }  
  void Print()const{
    printf("%2.3f %2.3f %2.3f %2.3f %2.3f\t", cmean, csig, fmean, fsig, ffrac);
    //    cout<<cmean<<"\t"<<csig<<"\t"
    //<<fmean<<"\t"<<fsig<<"\t"<<ffrac<<"\t";
  }

  void Write(ofstream &outfile)const{
    //printf("%2.3f %2.3f %2.3f %2.3f %2.3f\t", cmean, csig, fmean, fsig, ffrac);
    outfile<<cmean<<"\t"<<csig<<"\t"
       <<fmean<<"\t"<<fsig<<"\t"<<ffrac<<"\t";
  }
};
#undef GRND
#endif



int main(int argc, char**argv){
  cout<<r0ofmag(-19   )<<" "<<gamofmag(-19)<<" ";
  cosmo.Print();
  cout<<sim.Boxsize()<<endl;
  cout<<sim.Boxsize()*BOXFR<<endl;
  cosmo.ReadZFile();
  if(argc!=2){
    cout<<"Usage:pdfxi filelabel"<<endl;
    exit(1);
  }
  string label = argv[1];
  //  string cftmpfn = "./gal_cf"+label+".dat";
  string outfn = "./output_"+label+"_31k.dat";

  string com1 = "mkdir cfdat_"+label;
  system(com1.c_str());
  // Read in particles from cubes.
  vector <Particle*> particles = ReadParticles();
  cout<<" Read "<<particles.size()<<" particles"<<endl;
  assert(particles.size()>0);

  cout<<"Sorting"<<endl;
  // Sort particles on local mass density
  sort(particles.begin(),particles.end(),DLess);

  ofstream datfile(outfn.c_str());


  //gamma version
  //note::"sig" is actually the variance.

  int i = 0;
#ifdef LOOP
  cout<<"Looping"<<endl;
  for(double cmean=-0.3;cmean>=-0.71;cmean-=0.2){
    for(double csig=0.46;csig>=0.3399;csig-=0.06){
      for(double fmean=1.0;fmean<=3.501;fmean+=0.25){
	for(double fsig=0.4;fsig<=2.405;fsig+=0.25){
	  for(double ffrac=0.2;ffrac<=1.01;ffrac+=0.1){
#endif
#ifdef LIST
	    cout<<"Running list"<<endl;
	    string filename= "denspdf/hv_denspdf.dat";
	    ifstream file(filename.c_str());
	    if (file.fail()) {
	      cerr<<"error: cannot open "<<filename<<endl;
	      exit(1);
	    }
	    
	    int entries=28;//100;//00;//2690;
	    vector <double> cm(entries);
	    vector <double> cs(entries);
	    vector <double> fm(entries);
	    vector <double> fs(entries);
	    vector <double> ff(entries);
	    vector <double> mr(entries);
	    for (int ii=0; ii<entries; ii++) {
	      file>>mr[ii]>>cm[ii]>>fm[ii]>>fs[ii]>>ff[ii];
	    }
	    file.close();	      

	    for (int ii=0; ii<entries; ii++) {
	      cout<<mr[ii]<<" "<<r0ofmag(mr[ii])<<" "<<gamofmag(mr[ii])<<" ";
	      double cmean = cm[ii];
	      double csig = 0.4;//0.35;//cs[ii];
	      double fmean = fm[ii];
	      double fsig = fs[ii];
	      double ffrac = ff[ii];
	      cout<<cmean<<"\t"<<csig<<"\t"<<fmean<<"\t"<<fsig<<"\t"<<ffrac<<" ";	      
	      cout<<r0<<" "<<gam<<" "<<acctot<<" "<<chi2<<endl;
#endif	      
#ifdef MCMC
	     float acctot = 0.;
	      int Nchain = 10000;
	      //  string outmcmc = "mcmc_chain.dat";
	      double oldlik=0, newlik=0, accept=0;
	      ChainEl oldel, newel;
	      string chain_out = "mcmc_chain_chi_"+label+simlabel+".dat";
	      ofstream chout(chain_out.c_str());
	      for(int i=0; i<Nchain; i++){
		runcount++;
		//cout<<i<<" ";
#endif
		double r0, gam, chi2;

		pdf_to_xi(particles, newel.cmean, newel.csig, newel.fmean, newel.fsig, newel.ffrac, r0, gam, chi2, datfile, label);
		
#ifdef MCMC
		bool acc = false;
		accept = 0;
		if(i==0){
		  oldlik = -0.5*chi2;
		}
		else{
		  newlik = -0.5*chi2;
		  if (newlik > oldlik){
		    accept = 1.0;
		    acc = 1;
		    //cout<<"h "<<newlik<<" "<<oldlik<<" "<<acc<<endl;
		  }
		  else
		    accept = exp(newlik-oldlik);
		}
		if(drand48()<accept){
		  acc = 1;
		  //cout<<"h2 "<<accept<<" "<<acc<<endl;
		}
		//cout<<"h3 "<<accept<<" "<<acc<<endl;
		//if((gam<-2.1)||(r0<4.2)||(r0>12)||(gam>-1.7)) 
		// acc = false;
		//		chout<<" "<<oldlik<<" "<<newlik<<" \t"<<r0<<" "<<gam<<" "<<chi2<<" "<<accept<<" "<<acc<<endl;		      
		chout<<i<<" "<<acc<<" ";
		newel.Write(chout);
		chout<<chi2<<endl;
		acctot +=acc;
		cout<<i<<" "<<acc<<" ";
		printf("%7.4f %7.4f %7.4f %7.4f %7.4f\t", newel.cmean, newel.csig, newel.fmean, newel.fsig, newel.ffrac);
		//		cout<<newel.cmean<<"\t"<<newel.csig<<"\t"<<newel.fmean<<"\t"<<newel.fsig<<"\t"
		//  <<newel.ffrac<<"\t";	      
		cout<<r0<<" "<<gam<<" "<<acctot/(i+1.)<<" "<<oldlik<<" "<<newlik<<" "<<chi2<<endl;
		if (acc) {
		  oldel = newel;
		  oldlik = newlik;
		  //cout<<"resetting"<<accept<<" "<<oldel.fmean<<" ";
		  //cout<<newel.fmean<<endl;
		  //newel.Write(chout);
		}
		newel.Reset(oldel);
	      }
#endif
#ifdef WRITECF
	      string com = "cp "+cffile+" "+"cfdat_"+label+"/cf"+MakeString(i,4)+".dat";
	      system(com.c_str());
#endif
	      
	      i++;
#ifdef LOOP
	    }
	  }   	
	} 
      }
    }
#endif
#ifdef LIST
  }
#endif
  cout<<"done with program."<<endl;
}

