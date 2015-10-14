#ifndef sample_h
#define sample_h

enum sample_type{BRIGHT, MID, DIM, DIMMER, ZR20, ZR205, ZR21};
enum cf_type{XI, WP, WPCOVAR};
using namespace std;

class Sample{
 public:
  Sample(sample_type samplet, int tnbin, int tstart, int tend, cf_type cft):
    nbin(tnbin), startpoint(tstart), endpoint(tend), cftype(cft){
    if(samplet==BRIGHT){
      magmin=-25;
      magmax = -22;
      if(cft==XI)
	xi_filename = "xirk_23_22.dat";
      else
	xi_filename = "wp_23_22.dat";
      xilabel = "cfb";
      nbins = 9;
      rmin = 0.812831;
      rmax = 20.4174;
      ncf = 20000;
      //cout<<nbin<<" "<<ncf<<endl;
    }
    if(samplet==MID){
      magmin=-22;
      magmax = -21;
      if(cft==XI)
      xi_filename = "xirk_22_21.dat";
      else
	xi_filename = "wp_22_21.dat";
      xilabel = "cfm";
      ncf = 23930;
      nbins = 12;
      rmin = 0.128825;
      rmax = 20.4174;
    }
    if(samplet==DIM){
      magmin=-21;
      magmax = -20;
      if(cft==XI)
	xi_filename = "xirk_21_20.dat";
      else
	xi_filename = "wp_21_20.dat";
      xilabel = "cfd";
      ncf = 31053;
      nbins = 12;
      rmin = 0.128825;
      rmax = 20.4174;
    }
  }

  float Magmin()const{return magmin;};
  float Magmax()const{return magmax;};
  float Rmin()const{return rmin;};
  float Rmax()const{return rmax;};
  int Nbins()const{return nbins;};
  int Nbin()const{return nbin;};
  int Ncf()const{return ncf;};
  int Startpoint()const{return startpoint;};
  int Endpoint()const{return endpoint;};
  string Simfile(string dir, int runcout){
    string run = MakeString(runcout, 5);
    string s = dir+xilabel+run+".dat";
    return s;
  }
  string Simfile(string dir, int runcout, int iter){
    string run = MakeString(runcout, 5);
    string it = MakeString(iter, 1);
    string s = dir+xilabel+run+'.'+it+".dat";
    return s;
  }
  string Datfile(){
    return xi_filename;
  }
  cf_type Cftype(){
    return cftype;
  }
 private:
  int nbins;
  int nbin;
  int ncf;
  int startpoint;
  int endpoint;
  float magmin;
  float magmax;
  float rmin;
  float rmax;
  string xilabel;
  string xi_filename;
  cf_type cftype;
};

#endif
