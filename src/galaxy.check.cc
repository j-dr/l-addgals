#include <cmath>
#include <iostream>
#include "nr.h"
#include "myrand.h"
#include "fivetuple.h"
#include "constants.h"
#include "galaxy.h"
#include "cosmo.h" 
#ifdef LF_FROM_DATA
#include <vector>
#endif

using namespace std;

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

//#define NBIN 85000
#define NBIN 850

struct den_ent{
  float x[NBIN];
};

float LocalDens(FiveTuple fp, FiveTuple fpp, int weight1, int weight2)
{
  //generate a table of P(y) = x
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

float LocalDens(den_ent x, den_ent y)
{
  //for(int i=0;i<NBIN;i++)
  //cout<<x.x[i]<<" "<<y.x[i]<<endl;
  float d8_max = 8.5;
  float d8_min = 0.1;
  int nbin = NBIN;
  float d8 = d8_min - 1.0;
  while(d8 < d8_min || d8 > d8_max){
    float ranu = drand48();
    int ind = 0;
    for(ind=0;ind<nbin;ind++){
      if(ranu < x.x[ind])
	break;
    }
    if (ind == 0)
      ind++;
    if (ind == nbin)
      //nbin++;
      ind--;
    //cout<<" ranu = "<<ranu<<" ind = "<<ind<<" x[ind] = "<<x.x[ind];
    float dx = ranu - x.x[ind-1];
    float slope = (y.x[ind]-y.x[ind-1])/(x.x[ind]-x.x[ind-1]);
    d8 = y.x[ind-1]+dx*slope;
    //cout<<" d8 = "<<d8<<endl;
    //cout<<"Diagnistics: "<<endl;
    //cout<<"    "<<dx<<" "<<ranu<<" "<<x.x[ind-1]<<endl;
    //cout<<"    "<<slope<<" "<<y.x[ind]<<" "<<y.x[ind-1]<<" "<<x.x[ind]<<" "<<x.x[ind-1]<<endl;

  }
  if (d8 > d8_max)
    d8 = d8_max;
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
  float gmagbin1, gdensity1;
  ifstream file(lf_file.c_str());
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
 //return gdensity[i-1] + diff*slope;
 return LF;
  //return LumNumberDensityInterp(M);
#else
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
#ifdef MAG_LIMITED
  double mmax = -1.;
#else
  double mmax = -8.;
#endif
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

void make_denspdf(vector <FiveTuple> &denspdf, vector <float> &mags,ChainEl chel){
  float max_mag = -22.8;
  float min_mag = -18.0;
  float del_mag = 0.1;
  float tmag = min_mag;
  ofstream file_out("denspdf_used.ascii");
  while(tmag > max_mag){
    mags.push_back(tmag);
    float lum = pow(10.,-0.4*(tmag-Mstar));
    FiveTuple fTup(chel.cmean(),chel.fmean(lum),chel.fsig(lum),chel.ffrac(lum)); 
    file_out<<tmag<<" "<<fTup[0]<<" "<<fTup[1]<<" "<<fTup[2]<<" "<<fTup[3]<<" "<<fTup[4]<<endl;
    denspdf.push_back(fTup);
    tmag -= del_mag;
  }
}

void ReadPDFs(vector <FiveTuple> &denspdf, vector <float> &mags){
  string filename = "denspdf/"+simlabel+"_denspdf.dat";

#ifdef EVOLVE_DENSPDF
  //#ifndef SNAPSHOT
  /*
  if (ZREDMIN > 0.1)
    filename = "denspdf/"+simlabel+"_denspdf_088.dat";
  if (ZREDMIN > 0.5)
    filename = "denspdf/"+simlabel+"_denspdf_080.dat";
  if (ZREDMIN > 0.8)
    filename = "denspdf/"+simlabel+"_denspdf_072.dat";
  */

  //string pdf_base = "Consuelo02_all_noscatter_montero_denspdf";
  string pdf_base = "Consuelo02_centrals_scatter_montero_denspdf";
  filename = "denspdf/"+pdf_base+"_099.dat";
  if (ZREDMIN > 0.027)
    filename = "denspdf/"+pdf_base+"_098.dat";
  if (ZREDMIN > 0.054)
    filename = "denspdf/"+pdf_base+"_097.dat";
  if (ZREDMIN > 0.082)
    filename = "denspdf/"+pdf_base+"_096.dat";
  if (ZREDMIN > 0.11)
    filename = "denspdf/"+pdf_base+"_095.dat";
  if (ZREDMIN > 0.14)
    filename = "denspdf/"+pdf_base+"_094.dat";
  if (ZREDMIN > 0.17)
    filename = "denspdf/"+pdf_base+"_093.dat";
  if (ZREDMIN > 0.2)
    filename = "denspdf/"+pdf_base+"_092.dat";
  if (ZREDMIN > 0.23)
    filename = "denspdf/"+pdf_base+"_091.dat";
  if (ZREDMIN > 0.27)
    filename = "denspdf/"+pdf_base+"_090.dat";
  if (ZREDMIN > 0.30)
    filename = "denspdf/"+pdf_base+"_089.dat";
  if (ZREDMIN > 0.33)
    filename = "denspdf/"+pdf_base+"_088.dat";
  if (ZREDMIN > 0.37)
    filename = "denspdf/"+pdf_base+"_087.dat";
  if (ZREDMIN > 0.41)
    filename = "denspdf/"+pdf_base+"_086.dat";
  if (ZREDMIN > 0.44)
    filename = "denspdf/"+pdf_base+"_085.dat";
  if (ZREDMIN > 0.48)
    filename = "denspdf/"+pdf_base+"_084.dat";
  if (ZREDMIN > 0.52)
    filename = "denspdf/"+pdf_base+"_083.dat";
  //if (ZREDMIN > 0.56)
  //filename = "denspdf/"+pdf_base+"_082.dat";
  if (ZREDMIN > 0.60)
    filename = "denspdf/"+pdf_base+"_081.dat";
  if (ZREDMIN > 0.64)
    filename = "denspdf/"+pdf_base+"_080.dat";
  if (ZREDMIN > 0.69)
    filename = "denspdf/"+pdf_base+"_079.dat";
  if (ZREDMIN > 0.73)
    filename = "denspdf/"+pdf_base+"_078.dat";
  if (ZREDMIN > 0.78)
    filename = "denspdf/"+pdf_base+"_077.dat";
  if (ZREDMIN > 0.83)
    filename = "denspdf/"+pdf_base+"_076.dat";
  if (ZREDMIN > 0.87)
    filename = "denspdf/"+pdf_base+"_075.dat";
  if (ZREDMIN > 0.92)
    filename = "denspdf/"+pdf_base+"_074.dat";
  if (ZREDMIN > 0.97)
    filename = "denspdf/"+pdf_base+"_073.dat";
  if (ZREDMIN > 1.03)
    filename = "denspdf/"+pdf_base+"_072.dat";
  if (ZREDMIN > 1.08)
    filename = "denspdf/"+pdf_base+"_071.dat";
  if (ZREDMIN > 1.14)
    filename = "denspdf/"+pdf_base+"_070.dat";
  if (ZREDMIN > 1.19)
    filename = "denspdf/"+pdf_base+"_069.dat";
  if (ZREDMIN > 1.25)
    filename = "denspdf/"+pdf_base+"_068.dat";
  if (ZREDMIN > 1.31)
    filename = "denspdf/"+pdf_base+"_067.dat";
  if (ZREDMIN > 1.37)
    filename = "denspdf/"+pdf_base+"_066.dat";
  if (ZREDMIN > 1.50)
    filename = "denspdf/"+pdf_base+"_064.dat";
  if (ZREDMIN > 1.63)
    filename = "denspdf/"+pdf_base+"_062.dat";
  //#endif
#endif

#ifdef NON_PARAMETERIZED_PDF
  //Read the Dimensions
  file>>
  while(file)

  return
#endif

  ifstream file(filename.c_str());
  if (file.fail()) {
    cerr<<"error: cannot open "<<filename<<endl;
    exit(1);
  }
  else
    cerr<<"reading "<<filename<<endl;

  float mag, cm, cs, fm, fs, p;
  int nparams = 4;

  file>>mag>>cm>>cs>>fm>>fs>>p;
  if (p > 0){  //kludgy-way of seing if this is a 4 or 5 parameter PDF
    //cout<<"Reading file '"<<simlabel<<"_denspdf.dat' with magniutude + 5 parameter PDF"<<endl;
    cout<<"Reading file '"<<filename<<"' with magniutude + 5 parameter PDF"<<endl;
    nparams = 5;
    FiveTuple fTup(cm,cs,fm,fs,p);
    mags.push_back(mag);
    denspdf.push_back(fTup);
    fTup.Print();
    //cout<<mag<<" "<<cm<<" "<<cs<<" "<<fm<<" "<<fs<<" "<<p<<endl;
    while(file){
      file>>mag>>cm>>cs>>fm>>fs>>p;
      //cout<<mag<<" "<<cm<<" "<<cs<<" "<<fm<<" "<<fs<<" "<<p<<endl;
      FiveTuple fTup(cm,cs,fm,fs,p);
      fTup.Print();
      mags.push_back(mag);
      denspdf.push_back(fTup);
    }   
  } else{
    //cout<<"Reading file '"<<simlabel<<"_denspdf.dat' with magniutude + 4 parameter PDF"<<endl;
    cout<<"Reading file '"<<filename<<"' with magniutude + 4 parameter PDF"<<endl;
    file.clear();
    file.seekg(0,ios::beg);

    cout<<endl;

    while(file){
      //float tmp;
      //float mag, cm, fm, fs, p;
      //file>>mag>>tmp>>tmp>>cm>>fm>>fs>>p>>tmp>>tmp;
      file>>mag>>cm>>fm>>fs>>p;
      //cout<<mag<<" "<<cm<<" "<<fm<<" "<<fs<<" "<<p<<endl;
      FiveTuple fTup(cm,fm,fs,p);
      //fTup.Print();
      mags.push_back(mag);
      denspdf.push_back(fTup);
    }

    cout<<endl;

  }

  /*
  //This part just smoothes the fm, fs, p values
  for(int i=1;i<denspdf.size()-1;i++){
    float cm, cs, fm, fs, p;
    //    denspdf[i].Print();
    cm = denspdf[i][0];
    cs = denspdf[i][1];
    fm = (denspdf[i-1][2]+denspdf[i][2]+denspdf[i+1][2])/3.;
    fs = (denspdf[i-1][3]+denspdf[i][3]+denspdf[i+1][3])/3.;
    p = (denspdf[i-1][4]+denspdf[i][4]+denspdf[i+1][4])/3.; 
   //FiveTuple fTup(cm,fm,fs,p); 
    FiveTuple fTup(cm,cs,fm,fs,p); 
    denspdf[i] = fTup;
    //    denspdf[i].Print();
    //    for(int j=2;j<=4;j++){
      //      cout<<j<<" "<<denspdf[i-1][j]<<" "
      //  <<denspdf[i][j]<<" "
      //  <<denspdf[i+1][j]<<endl;
      //float newval = denspdf[i][j];
      
      //}
  }
  */

  cout<<"Read "<<denspdf.size()<<"lines from denspdf"<<endl;;
}

void define_prob(struct den_ent * &x_prob, struct den_ent * &y_prob, vector <int> &weight, vector <FiveTuple> &denspdf, vector <float> &dmagbins, float vol)
{
  for(int i=0;i<dmagbins.size();i++){
    int tweight = (int) (LumNumberDensity(dmagbins[i])*vol);
    weight.push_back(tweight);
  }  
  //x_prob = new den_ent [NBIN];
  //y_prob = new den_ent [NBIN];
  int nbin = NBIN;
  float y_del = 8.5/nbin;
  for(int id=0;id<denspdf.size()-1;id++){
    den_ent x;
    den_ent y;
    float weight1 = weight[id];
    float weight2 = weight[id+1];
    FiveTuple fp = denspdf[id];
    FiveTuple fpp = denspdf[id+1];
    //fp.Print();
    for(int i=0;i<nbin;i++){
      y.x[i] = (i+1)*y_del;
      float p1 = 0.5*(1. - fp[4])*(1+erf((log(y.x[i])-fp[0])/(fp[1]*sqrt(2.0))));
      float p2 = 0.5*fp[4]*(1+erf((y.x[i]-fp[2])/(fp[3]*sqrt(2.0))));
      float p3 = 0.5*(1. - fpp[4])*(1+erf((log(y.x[i])-fpp[0])/(fpp[1]*sqrt(2.0))));
      float p4 = 0.5*fpp[4]*(1+erf((y.x[i]-fpp[2])/(fpp[3]*sqrt(2.0))));
      x.x[i] = weight1*(p1+p2) - weight2*(p3+p4);
    }
    //float offset = x.x[0];
    //for(int i=0;i<nbin;i++)
    //x.x[i] -= offset;
    float norm = x.x[nbin-1];
    if (norm > 0){
      for(int i=0;i<nbin;i++)
	x.x[i] /= norm;
    }
    x_prob[id] = x;
    y_prob[id] = y;
  }
  ofstream pdf_file("pdf_test.ascii");
  for(int id=0;id<nbin;id++){
    pdf_file<<x_prob[36].x[id]<<" "<<y_prob[36].x[id]<<endl;
  }
  //cout<<"Some x 0's: ";
  //for(int i=0;i<denspdf.size();i++){
  //cout<<" "<<x_prob[i].x[0]<<" ";
  //}
  cout<<endl;
}

/*
vector <Galaxy *> GetGalaxies(double vol, 
		 float ff, float cm, float cs){
  integrator LF_integrator(&LF, 0, 1.0E-05, 0.01);
  int n = (int) (LumNumberDensity(Magmin)*vol);
  vector <Galaxy*> galaxies;
  galaxies.reserve(n);
  double volume = vol;
  int ngal = 0;

  for(double mag=-23.7; mag <=-16; mag+=0.001){
      double num = LF_integrator.Integrate(-30, mag)*volume;
      int diff = ((int) floor(num))-galaxies.size();
      if (num>n) diff = n-galaxies.size();
      for(int i=0;i<diff;i++){
	Galaxy * galaxy = new Galaxy(mag,ngal, ff, cm, cs);//, fms, fmi, fss, fsi);
 	galaxies.push_back(galaxy);
	ngal++;
      }
      if (galaxies.size()>=n) break;
  }  
  return galaxies;
}
*/
vector <Galaxy *> GetGalaxies(double vol){
  integrator LF_integrator(&LF, 0, 1.0E-05, 0.01);
  cout<<vol<<endl;
  int n = (int) (LumNumberDensity(Magmin)*vol);
  PRNT("GetGalaxies:", n);
  vector <Galaxy*> galaxies;
  galaxies.reserve(n);
  double volume = vol;
  cout<<volume<<" "<<sim.CubeVolume()<<endl;
  int ngal = 0;
  ofstream dens_file("dens_check.ascii");
  ofstream denspdf_file("denspdf_check.ascii");
  ofstream gal_prop("galaxies.ascii");

  vector <FiveTuple> denspdf;
  vector <float> dmagbins;
  ReadPDFs(denspdf, dmagbins); 
  assert(denspdf.size()== dmagbins.size());
  float pdfbinsize = dmagbins[0]-dmagbins[1];
  int ind = dmagbins.size()-3;
  //int ind = dmagbins.size()-4;
  FiveTuple fTup = denspdf[ind];
  FiveTuple fTup_prime = denspdf[ind+1];
  int weight1, weight2;
  float cur_mag = dmagbins[ind];
  cout<<"Checking the first denspdf values: ";
  denspdf[0].Print();
  //We need to calculate the total number of particles for each bin
  /*
  vector <int> weight;
  for(int i=0;i<dmagbins.size();i++){
    int tweight = (int) (LumNumberDensity(dmagbins[i])*vol);
    weight.push_back(tweight);
  } 
  */ 
  cout<<" The pdf parameters (has "<<dmagbins.size()<<" elements)"<<endl;

  //for(int i=0;i<dmagbins.size();i++){
    //cout<<weight[i]<<" "<<dmagbins[i]<<" ";
    //denspdf[i].Print();
  //}

  cout<<endl;

#ifdef DENSPDF_IN_MAG_BIN
  struct den_ent * x_prob;
  struct den_ent * y_prob;
  x_prob = new den_ent [NBIN];
  y_prob = new den_ent [NBIN];
  vector <int> weight;
  cout<<"Defining the probability arrays."<<endl;
  define_prob(x_prob, y_prob, weight, denspdf, dmagbins, vol);
  cout<<"done."<<endl;
  fTup = denspdf[ind];
  fTup_prime = denspdf[ind+1];
  weight1 = weight[ind];
  weight2 = weight[ind+1];
  cur_mag = dmagbins[ind];
  for(int i=0;i<dmagbins.size();i++)
    denspdf_file<<weight[i]<<" "<<dmagbins[i]<<endl;
  /*
  struct den_ent * x_prob;
  x_prob = new den_ent [NBIN];
  struct den_ent * y_prob;
  y_prob = new den_ent [NBIN];
  int nbin = NBIN;
  float y_del = 8.5/nbin;
  for(int id=0;id<dmagbins.size()-1;id++){
    den_ent x;
    den_ent y;
    float weight1 = weight[id];
    float weight2 = weight[id+1];
    FiveTuple fp = denspdf[id];
    FiveTuple fpp = denspdf[id+1];
    for(int i=0;i<nbin;i++){
      y.x[i] = (i+1)*y_del;
      float p1 = 0.5*(1. - fp[4])*(1+erf((log(y.x[i])-fp[0])/(fp[1]*sqrt(2.0))));
      float p2 = 0.5*fp[4]*(1+erf((y.x[i]-fp[2])/(fp[3]*sqrt(2.0))));
      float p3 = 0.5*(1. - fpp[4])*(1+erf((log(y.x[i])-fpp[0])/(fpp[1]*sqrt(2.0))));
      float p4 = 0.5*fpp[4]*(1+erf((y.x[i]-fpp[2])/(fpp[3]*sqrt(2.0))));
      x.x[i] = weight1*(p1+p2) - weight2*(p3+p4);
    }
    //float offset = x.x[0];
    //for(int i=0;i<nbin;i++)
    //x.x[i] -= offset;
    float norm = x.x[nbin-1];
    for(int i=0;i<nbin;i++)
      x.x[i] /= norm;
    x_prob[id] = x;
    y_prob[id] = y;
  }
  */
  //ofstream pdf_file("pdf_test.ascii");
  //for(int id=0;id<NBIN;id++){
  //pdf_file<<x_prob[36].x[id]<<" "<<y_prob[36].x[id]<<endl;
  //}
#endif

  // cout<<cur_mag<<" "<<ind<<" "<<dmagbins[ind-1]<<" "<<dmagbins[ind-2]<<endl;
  /*
  ofstream testout("testdens.dat");
  fTup = denspdf[1];
  fTup.Print();
  for(int t=0;t<1000;t++){
    testout<<fTup.LocalDens()<<endl;
  }

  ofstream testout1("testdens1.dat");
  FiveTuple fTup1(-0.63,3.48,0.96,1.0);
  fTup1.Print();
  for(int t=0;t<1000;t++){
    testout1<<fTup1.LocalDens()<<endl;
  }

  ofstream testout2("testdens2.dat");

  FiveTuple fTup2(-0.63,3.48,0.96,0.0);
  fTup2.Print();
  for(int t=0;t<1000;t++){
    testout2<<fTup2.LocalDens()<<endl;
  }
  */

  cout<<"Generating Galaxy Magnitudes."<<endl;
  double ng_expected = LF_integrator.Integrate(-30, Magmin)*volume;
  cout<<"Expect to generate "<<ng_expected<<" galaxies."<<endl;
  for(double mag=-23.7; mag <=-8; mag+=0.001){
    //cout<<"mag = "<<mag<<", cur_mag = "<<cur_mag<<endl;
    //if(mag>cur_mag+0.5*pdfbinsize){
    if(mag>cur_mag){
      if(ind >= 0){
	if (ind != 0)//mbusha added if clause to get proper dim binning
	  --ind;
	fTup = denspdf[ind];
	fTup_prime = denspdf[ind+1];
	weight1 = weight[ind];
	weight2 = weight[ind+1];
	cur_mag = dmagbins[ind];
	//cout<<ind<<" "<<mag<<" "<<cur_mag<<" "<<cur_mag+0.5*pdfbinsize<<" "<<galaxies.size()<<" "<<n<<endl;
      }
      else{
	cout<<"[GetGalaxies] PDF binning problem... exiting"<<endl;
	cout<<ind<<" "<<mag<<" "<<cur_mag<<" "<<cur_mag+0.5*pdfbinsize<<endl;
	exit(1);
      }
    }
    double num = LF_integrator.Integrate(-30, mag)*volume;
    int diff = ((int) floor(num))-galaxies.size();
    if (num>n) diff = n-galaxies.size();
    //cout<<"Using ind = "<<ind<<", going to add "<<diff<<" new galaxies.  Currently have "<<galaxies.size()<<" galaxiees."<<endl;
    for(int i=0;i<diff;i++){
      //cout<<"Getting galaxy redshift."<<endl;
      float zGal = SelectGalaxyZ();
      //cout<<"Selected z = "<<zGal<<endl;
#ifdef SNAPSHOT
      float ThisMag = mag;
      float ThisMStar = Mstar;
#else
      float ThisMag = evolve_mag(mag, zGal);
      float ThisMStar = evolve_mag(Mstar, zGal);
#endif
      float lum = pow(10.,-0.4*(ThisMag-ThisMStar));
      //cout<<"have lum = "<<lum<<endl;
#ifdef DENSPDF_IN_MAG_BIN
      //float dist8 = LocalDens(fTup, fTup_prime, weight1, weight2);
      //cout<<"Calculating dist 8 for ind = "<<ind<<endl;
      //for(int id=0;id<NBIN;id++){
      //cout<<x_prob[ind].x[id]<<" "<<y_prob[ind].x[id]<<endl;
      //}
      //cout<<"  Getting LocalDens"<<endl;
      float dist8 = LocalDens(x_prob[ind], y_prob[ind]);
      //cout<<"Got dist8 = "<<dist8<<endl;
      dens_file<<mag<<" "<<zGal<<" "<<dist8<<" "<<ind<<endl;
      Galaxy * galaxy = new Galaxy(mag,ngal,dist8);
#else
      Galaxy * galaxy = new Galaxy(mag,ngal,fTup.LocalDens());
#endif
      galaxies.push_back(galaxy);
      galaxies[ngal]->zGal(zGal);
      ngal++;
      gal_prop<<mag<<" "<<zGal<<" "<<dist8<<" "<<endl;
    }
    if (galaxies.size()>=n) break;
  }  
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


#ifdef MAG_LIMITED_chainel
vector <Galaxy *> GetDimGalaxies(double vol){
  integrator LF_integrator(&LF, 0, 1.0E-05, 0.01);
  int n = 0;
  double dz = 0.01;
  double zNow = ZREDMIN;
  while(zNow < ZREDMAX){
    double r1 = cosmo.RofZ(zNow);
    double r2 = cosmo.RofZ(zNow+dz);
    double vol = 4.0/3.0*PI*(r2*r2*r2-r1*r1*r1);
    double dl = (r1+0.0001)*1e6*(1+cosmo.ZofR(r1))*(1+cosmo.ZofR(r1));
    double magmin_z = oMagMin - 5*(log10(dl)-1);
    //cout<<" z = "<<zNow<<", magmin = "<<magmin_z<<", vol = "<<vol;
    int nnew = 0;
    if(magmin_z > Magmin)
      {
        //cout<<", adding galaxies....";
        nnew = (int) (LumNumberDensity(Magmin_pdf, magmin_z)*vol);
        n += nnew;
      }
    //cout<<", new galaxies = "<<nnew<<endl;
    zNow += dz;
  }

  PRNT("GetDimGalaxies:", n);
  vector <Galaxy*> galaxies;
  galaxies.reserve(n);
  double volume = vol;
  cout<<volume<<" "<<sim.CubeVolume()<<endl;
  int ngal = 0;
  double LowestMag;

  vector <FiveTuple> denspdf;
  vector <float> dmagbins;
  ReadPDFs(denspdf, dmagbins);
  assert(denspdf.size()== dmagbins.size());
  float pdfbinsize = dmagbins[0]-dmagbins[1];
  int ind = dmagbins.size()-2;
  FiveTuple fTup = denspdf[ind];
  float cur_mag = dmagbins[ind];

  //Second the Magnitude-limited part
  cout<<"Dim sample has "<<galaxies.size()<<" galaxies."<<endl;
  cout<<"Generating the Magnitude Limited Sample"<<endl;
  for(double mag=Magmin; mag <0; mag+=0.001){
    float DlMax = pow(10.0, 0.2*(oMagmin-mag)+1);
    DlMax *= 1.1; //More padding
    float zmax = 0.1;
    float rmax = cosmo.RofZ(zmax);
    float Dl = rmax*(1+zmax)*(1+zmax);
    while(Dl<DlMax){
      zmax += 0.1;
      rmax = cosmo.RofZ(zmax);
      if (zmax >= ZREDMAX){
        zmax = ZREDMAX;
        rmax = cosmo.RofZ(zmax);
        break;
      }
      Dl = rmax*1e6*(1+zmax)*(1+zmax);
    }
    zmax *= 1.1; //this factor is padding
    rmax = cosmo.RofZ(zmax);
    float PI = 3.1415926535897932384626433832795;
    float decmin_rad = PI*(90 - DECMAX)/180;
    float decmax_rad = PI*(90 - DECMIN)/180;
    float ramin_rad = PI*RAMIN/180;
    float ramax_rad = PI*RAMAX/180;
    float dra = ramax_rad - ramin_rad;
    float r_zmin = cosmo.RofZ(ZREDMIN);
    float volume_dim = (ramax_rad-ramin_rad)*(cos(decmin_rad)-cos(decmax_rad))*(pow(rmax,3)-pow(r_zmin,3))/3.;
    double diff = LF_integrator.Integrate(mag-0.001, mag)*volume_dim;
    for(int i=0;i<diff;i++){
      float zGal = SelectGalaxyZ(rmax);
#ifdef SNAPSHOT
      float ThisMag = mag;
      float ThisMStar = Mstar;
#else
      float ThisMag = evolve_mag(mag, zGal);
      float ThisMStar = evolve_mag(Mstar, zGal);
#endif
      float dl = cosmo.RofZ(zGal)*1e6*(1+zGal)*(1+zGal);
      float omag = ThisMag + 5*(log10(dl)-1);  //causes a mismatch between bright and dim. Why?
      //float omag = mag + 5*(log10(dl)-1);
      //if (omag > oMagMin && mag > Magmin) //Why did I do this before???
      if (omag > oMagMin || mag < Magmin)
        continue;

      float lum = pow(10.,-0.4*(ThisMag-ThisMStar));
      FiveTuple fTup(chel.cmean(),chel.fmean(lum),chel.fsig(lum),chel.ffrac(lum));
      Galaxy * galaxy = new Galaxy(mag,ngal,fTup.LocalDens());
      galaxy->zGal(zGal);
      galaxies.push_back(galaxy);
      ngal++;
    }
    LowestMag = mag;
    if (galaxies.size()>=n) break;
  }

  cout<<"[GetDimGalaxies] Done generating "<<galaxies.size()<<" down to magnitude "<<LowestMag<<endl;
  return galaxies;
}
#endif

vector <Galaxy *> GetGalaxiesNodens(double vol){
  integrator LF_integrator(&LF, 0, 1.0E-05, 0.01);
  cout<<vol<<endl;
  int n = (int) (LumNumberDensity(Magmin_pdf)*vol);
  PRNT("GetGalaxies:", n);
  vector <Galaxy*> galaxies;
  galaxies.reserve(n);
  double volume = vol;
  cout<<volume<<" "<<sim.CubeVolume()<<endl;
  int ngal = 0;
                                                        
  for(double mag=-23.7; mag <=-8; mag+=0.001){

      double num = LF_integrator.Integrate(-30, mag)*volume;
      int diff = ((int) floor(num))-galaxies.size();
      if (num>n) diff = n-galaxies.size();
      for(int i=0;i<diff;i++){
	Galaxy * galaxy = new Galaxy(mag,ngal);
 	galaxies.push_back(galaxy);
	ngal++;
      }
      if (galaxies.size()>=n) break;
  }  
  cout<<"[GetGalaxies] Done assigning galaxies "<<galaxies.size()<<endl;
  return galaxies;
}

vector <Galaxy *> GetGalaxies(double vol, ChainEl chel){
  integrator LF_integrator(&LF, 0, 1.0E-05, 0.01);
  cout<<vol<<endl;

#ifdef DENSPDF_IN_MAG_BIN
  vector <FiveTuple> denspdf;
  vector <float> dmagbins;
  cout<<"Galling make_denspdf"<<endl;
  make_denspdf(denspdf, dmagbins, chel);
  cout<<"          done"<<endl;
  struct den_ent * x_prob;
  struct den_ent * y_prob;
  x_prob = new den_ent [NBIN];
  y_prob = new den_ent [NBIN];
  vector <int> weight;
  cout<<"Defining the probability arrays."<<endl;
  define_prob(x_prob, y_prob, weight, denspdf, dmagbins, vol);
  cout<<"done."<<endl;
  //for(int i=0;i<dmagbins.size();i++)
  //cout<<dmagbins[i]<<endl;
  ofstream d8_file("dist8_check.ascii");
#endif

  int n = (int) (LumNumberDensity(Magmin_pdf)*vol);//the 1.2 is padding
  PRNT("GetGalaxies:", n);
  vector <Galaxy*> galaxies;
  galaxies.reserve(n);
  double volume = vol;
  cout<<volume<<" "<<sim.CubeVolume()<<endl;
  int ngal = 0;
  double LowestMag;
  for(double mag=-23.7; mag <=-8; mag+=0.001){
      double num = LF_integrator.Integrate(-30, mag)*volume;
      int diff = ((int) floor(num))-galaxies.size();
      if (num>n) diff = n-galaxies.size();
      //float lum = pow(10.,-0.4*(mag+20.44));
      //FiveTuple fTup(chel.cmean(),chel.fmean(lum),chel.fsig(lum),chel.ffrac(lum)); 
      for(int i=0;i<diff;i++){
	float zGal = SelectGalaxyZ();
#ifdef SNAPSHOT
	float ThisMag = mag;
	float ThisMStar = Mstar;
#else
	float ThisMag = evolve_mag(mag, zGal);
	float ThisMStar = evolve_mag(Mstar, zGal);
#endif
	float lum = pow(10.,-0.4*(ThisMag-ThisMStar));
#ifdef DENSPDF_IN_MAG_BIN
	int ind = dmagbins.size()-2;
	while(dmagbins[ind] < ThisMag && ind > 0)
	  ind--;
	//cout<<"Looking for matnigude "<<ThisMag<<" with ind = "<<ind<<" total bins = "<<dmagbins.size()<<endl;
	//cout<<"  First few x_prob's: "<<x_prob[ind].x[0]<<" "<<x_prob[ind].x[1]<<" "<<x_prob[ind].x[2]<<" "<<x_prob[ind].x[3]<<" "<<x_prob[ind].x[4]<<endl;
	float dist8 = LocalDens(x_prob[ind], y_prob[ind]);
	Galaxy * galaxy = new Galaxy(mag,ngal,dist8);
	d8_file<<mag<<" "<<dist8<<endl;
#else
	FiveTuple fTup(chel.cmean(),chel.fmean(lum),chel.fsig(lum),chel.ffrac(lum)); 
	Galaxy * galaxy = new Galaxy(mag,ngal,fTup.LocalDens());
#endif
 	galaxies.push_back(galaxy);
	galaxies[ngal]->zGal(zGal);
	ngal++;
      }
      LowestMag = mag;
      if (galaxies.size()>=n) break;
  }  

  //  cout<<"[GetGalaxies] Done assigning galaxies "<<galaxies.size()<<endl;
  cout<<"[GetGalaxies] Done generating "<<galaxies.size()<<" down to magnitude "<<LowestMag<<endl;
  cout<<"[GetGalaxies] Range of magnitudes generated: "<<galaxies[0]->Mr()<<" "<<galaxies[galaxies.size()-1]->Mr()<<endl;
  return galaxies;
}

#ifdef MAG_LIMITED
vector <Galaxy *> GetDimGalaxies(double vol, ChainEl chel){
  integrator LF_integrator(&LF, 0, 1.0E-05, 0.01);

#ifdef DENSPDF_IN_MAG_BIN
  vector <FiveTuple> denspdf;
  vector <float> dmagbins;
  cout<<"Galling make_denspdf"<<endl;
  make_denspdf(denspdf, dmagbins, chel);
  cout<<"          done"<<endl;
  struct den_ent * x_prob;
  struct den_ent * y_prob;
  x_prob = new den_ent [NBIN];
  y_prob = new den_ent [NBIN];
  vector <int> weight;
  cout<<"Defining the dim probability arrays."<<endl;
  define_prob(x_prob, y_prob, weight, denspdf, dmagbins, vol);
  cout<<"done."<<endl;
#endif

  int n = 0;
  double dz = 0.01;
  double zNow = ZREDMIN;
  while(zNow < ZREDMAX){
    double r1 = cosmo.RofZ(zNow);
    double r2 = cosmo.RofZ(zNow+dz);
    double vol = 4.0/3.0*PI*(r2*r2*r2-r1*r1*r1);
    double dl = (r1+0.0001)*1e6*(1+cosmo.ZofR(r1))*(1+cosmo.ZofR(r1));
    double magmin_z = oMagMin - 5*(log10(dl)-1);
    //cout<<" z = "<<zNow<<", magmin = "<<magmin_z<<", vol = "<<vol;
    int nnew = 0;
    if(magmin_z > Magmin)
      {
	//cout<<", adding galaxies....";
	nnew = (int) (LumNumberDensity(Magmin_pdf, magmin_z)*vol);
	n += nnew;
      }
    //cout<<", new galaxies = "<<nnew<<endl;
    zNow += dz;
  }

  PRNT("GetDimGalaxies:", n);
  vector <Galaxy*> galaxies;
  galaxies.reserve(n);
  double volume = vol;
  cout<<volume<<" "<<sim.CubeVolume()<<endl;
  int ngal = 0;
  double LowestMag;
  
  cout<<"Dim sample has "<<galaxies.size()<<" galaxies."<<endl;
  cout<<"Generating the Magnitude Limited Sample"<<endl;
  for(double mag=Magmin; mag <0; mag+=0.001){
    float DlMax = pow(10.0, 0.2*(oMagmin-mag)+1);
    DlMax *= 1.1; //More padding
    float zmax = 0.1;
    float rmax = cosmo.RofZ(zmax);
    float Dl = rmax*(1+zmax)*(1+zmax);
    while(Dl<DlMax){
      zmax += 0.1;
      rmax = cosmo.RofZ(zmax);
      if (zmax >= ZREDMAX){
        zmax = ZREDMAX;
        rmax = cosmo.RofZ(zmax);
        break;
      }  
      Dl = rmax*1e6*(1+zmax)*(1+zmax);
    }
    zmax *= 1.1; //this factor is padding
    rmax = cosmo.RofZ(zmax);
    float PI = 3.1415926535897932384626433832795;
    float decmin_rad = PI*(90 - DECMAX)/180;
    float decmax_rad = PI*(90 - DECMIN)/180;
    float ramin_rad = PI*RAMIN/180;
    float ramax_rad = PI*RAMAX/180;
    float dra = ramax_rad - ramin_rad;
    float r_zmin = cosmo.RofZ(ZREDMIN);
    float volume_dim = (ramax_rad-ramin_rad)*(cos(decmin_rad)-cos(decmax_rad))*(pow(rmax,3)-pow(r_zmin,3))/3.;
    double diff = LF_integrator.Integrate(mag-0.001, mag)*volume_dim;
    for(int i=0;i<diff;i++){
      float zGal = SelectGalaxyZ(rmax);
      if (zGal > ZREDMAX || zGal < ZREDMIN)
	continue;
#ifdef SNAPSHOT
      float ThisMag = mag;
      float ThisMStar = Mstar;
#else
      float ThisMag = evolve_mag(mag, zGal);
      float ThisMStar = evolve_mag(Mstar, zGal);
#endif
      float dl = cosmo.RofZ(zGal)*1e6*(1+zGal)*(1+zGal);
      float omag = ThisMag + 5*(log10(dl)-1);  //causes a mismatch between bright and dim. Why?
      //float omag = mag + 5*(log10(dl)-1);
      //if (omag > oMagMin && mag > Magmin) //Why did I do this before???

      //I might want to remove this continue statement?
      if (omag > oMagMin || mag < Magmin)
	continue;
#ifdef DENSPDF_IN_MAG_BIN
      int ind = dmagbins.size()-2;
      while(dmagbins[ind] < ThisMag && ind > 0)
	ind--;
      float dist8 = LocalDens(x_prob[ind], y_prob[ind]);
      Galaxy * galaxy = new Galaxy(mag,ngal,dist8);
#else
      float lum = pow(10.,-0.4*(ThisMag-ThisMStar));
      FiveTuple fTup(chel.cmean(),chel.fmean(lum),chel.fsig(lum),chel.ffrac(lum)); 
      Galaxy * galaxy = new Galaxy(mag,ngal,fTup.LocalDens());
#endif
      galaxy->zGal(zGal);
      galaxies.push_back(galaxy);
      ngal++;
    }
    LowestMag = mag;
    if (galaxies.size()>=n) break;
  }  

  cout<<"[GetDimGalaxies] Done generating "<<galaxies.size()<<" down to magnitude "<<LowestMag<<endl;
  return galaxies;
}
#endif

/*
  this is wrong.
float area_coverage(float rmin, float rmax, float dmin, float dmax){
  float area = 0.;
  float interval = 0.001;
  for(float f = rmin; f<=rmax; f=f+interval){
    float sarea = (dmax-dmin)*interval;
    area = area+sarea*cos((f+0.5-90.)*3.141592654/180.);
  }
  return area;
}
*/

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


void GetDimGalaxies(vector <Galaxy*> &galaxies, double zmin, double zmax, double volume_frac, ChainEl chel){
  int nstart = 0;
  double vol = (4./3*PI*pow(cosmo.RofZ(zmax), 3)-4./3*PI*pow(cosmo.RofZ(zmin), 3))*fractional_area();
  //this is the smallest
  double amagmin = oMagmin-ZdistmodInterp(zmin);//-Q*(zmax-0.1);
  //this is the biggest, want it to be brighter for higher z
  double eMagmin = Magmin;//-Q*(zmax-0.1);
  if(amagmin>eMagmin){
    integrator LF_integrator(&LF, 0, 1.0E-05, 0.01);
    cout<<zmin<<" "<<zmin<<" "<<ZdistmodInterp(zmax)<<" "<<amagmin<<" "<<vol<<" "<<eMagmin<<endl;
    int n = (int) (LumNumberDensity(amagmin)*vol - LumNumberDensity(eMagmin));
    PRNT("GetGalaxies:", n);
    galaxies.reserve(n);
    double volume = vol;
    //  cout<<volume<<" "<<sim.CubeVolume()<<endl;
    int ngal = 0;
    
    for(double mag=eMagmin; mag <=0; mag+=0.001){
      double num = LF_integrator.Integrate(eMagmin, mag)*volume;
      int diff = ((int) floor(num))-(galaxies.size()-nstart);
      if (num>n) diff = n-(galaxies.size()-nstart);
      // evolve it back to z=0.1
      float lum = pow(10.,-0.4*(mag+20.44));//+(Q*zmax-0.1);
      FiveTuple fTup(chel.cmean(),chel.fmean(lum),chel.fsig(lum),chel.ffrac(lum)); 
      for(int i=0;i<diff;i++){
	Galaxy * galaxy = new Galaxy(mag,ngal,fTup.LocalDens());
 	galaxies.push_back(galaxy);
	ngal++;
      }
      if (galaxies.size()-nstart>=n) break;
    }  
  }
  cout<<"[GetGalaxies] Done assigning galaxies "<<galaxies.size()<<endl;
  return;
}



void define_prob_BCG(vector <Galaxy *> &galaxies, struct den_ent * &x_prob, struct den_ent * &y_prob, vector <int> &weight, vector <FiveTuple> &denspdf, vector <float> &dmagbins)
{
  for(int i=0;i<dmagbins.size();i++)
    weight.push_back(0);

  float nskip = 0;
  for(int i=0;i<galaxies.size();i++)
    {
      if (galaxies[i]->Central()){
	nskip++;
	continue;
      }
      for(int j=0;j<dmagbins.size();j++)
	if(galaxies[i]->Mr() < dmagbins[j] && galaxies[i]->Central() == 0)
	  weight[j]++;
    }
  cout<<"Skipped "<<nskip<<" centrals when calculating weight."<<endl;
	
  int nbin = NBIN;
  float y_del = 8.5/nbin;
  for(int id=0;id<denspdf.size()-1;id++){
    den_ent x;
    den_ent y;
    float weight1 = weight[id];
    float weight2 = weight[id+1];
    FiveTuple fp = denspdf[id];
    FiveTuple fpp = denspdf[id+1];
    for(int i=0;i<nbin;i++){
      y.x[i] = (i+1)*y_del;
      float p1 = 0.5*(1. - fp[4])*(1+erf((log(y.x[i])-fp[0])/(fp[1]*sqrt(2.0))));
      float p2 = 0.5*fp[4]*(1+erf((y.x[i]-fp[2])/(fp[3]*sqrt(2.0))));
      float p3 = 0.5*(1. - fpp[4])*(1+erf((log(y.x[i])-fpp[0])/(fpp[1]*sqrt(2.0))));
      float p4 = 0.5*fpp[4]*(1+erf((y.x[i]-fpp[2])/(fpp[3]*sqrt(2.0))));
      x.x[i] = weight1*(p1+p2) - weight2*(p3+p4);
    }
    float norm = x.x[nbin-1];
    if (norm > 0){
      for(int i=0;i<nbin;i++)
	x.x[i] /= norm;
    }
    x_prob[id] = x;
    y_prob[id] = y;
  }
}

void recalculate_d8(vector <Galaxy *> &galaxies){

  vector <FiveTuple> denspdf;
  vector <float> dmagbins;
  ReadPDFs(denspdf, dmagbins); 
  assert(denspdf.size()== dmagbins.size());
  float pdfbinsize = dmagbins[0]-dmagbins[1];
  int ind = dmagbins.size()-3;
  FiveTuple fTup = denspdf[ind];
  FiveTuple fTup_prime = denspdf[ind+1];
  int weight1, weight2;
  float mag;
  float cur_mag = dmagbins[ind];
  ofstream dens_file("redens_check.ascii");
  ofstream denspdf_file("redenspdf_check.ascii");

  struct den_ent * x_prob;
  struct den_ent * y_prob;
  x_prob = new den_ent [NBIN];
  y_prob = new den_ent [NBIN];
  vector <int> weight;
  define_prob_BCG(galaxies, x_prob, y_prob, weight, denspdf, dmagbins);
  fTup = denspdf[ind];
  fTup_prime = denspdf[ind+1];
  weight1 = weight[ind];
  weight2 = weight[ind+1];
  cur_mag = dmagbins[ind];

  for(int i=0;i<dmagbins.size();i++)
    denspdf_file<<weight[i]<<" "<<dmagbins[i]<<" "<<endl;

  for(int i=0;i<galaxies.size();i++){
    if(galaxies[i]->Central())
      continue;
    mag = galaxies[i]->Mr();
    
    //cout<<"Looking at galaxy "<<i<<" of "<<galaxies.size()<<" with mag = "<<mag<<endl;

    ind = dmagbins.size()-3;
    cur_mag = dmagbins[ind];
    weight1 = 0;
    weight2 = 0;
    //while(mag>cur_mag){
    while(mag>cur_mag | weight1 == weight2){
      if(ind > 0){   //if we're too dim we just stay in the last bin
	--ind;
	fTup = denspdf[ind];
	fTup_prime = denspdf[ind+1];
	weight1 = weight[ind];
	weight2 = weight[ind+1];
	cur_mag = dmagbins[ind];
      }
      else{
	break;
	/*
	cout<<"[GetGalaxies] PDF binning problem... exiting"<<endl;
	cout<<ind<<" "<<mag<<" "<<cur_mag<<" "<<cur_mag+0.5*pdfbinsize<<endl;
	exit(1);
	*/

      }
    }

    float dist8 = LocalDens(x_prob[ind], y_prob[ind]);
    galaxies[i]->Dist8(dist8);
    dens_file<<mag<<" "<<dist8<<" "<<ind<<" "<<galaxies[i]->Dist8()<<endl;
  }
}
