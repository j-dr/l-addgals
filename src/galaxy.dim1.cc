#include <cmath>
#include <iostream>
#include "nr.h"
#include "myrand.h"
#include "fivetuple.h"
#include "constants.h"
#include "galaxy.h"
#include "cosmo.h" 

using namespace std;

float FiveTuple::LocalDens() const{

  // returns a point selected from the
  // sum of a gaussian and a lognormal
  // with parameters specifed by the ftuple
  float min = -2.30259; //log(0.1)
  float d8;
  if(randbool(p)){
    float f8=-1;
    if((fm<0.1)||(fm>8.5)) cout<<"extreme f8"<<fm<<" "<<fs<<endl;
    while((f8<0.1)||(f8>8.5))
      f8 = fm+sqrt(fs)*NR::gasdev(seed);
    d8 = f8;
  }
  else{
    float c8 = min;
    if(cm<min) cout<<"extreme c8"<<cm<<" "<<cs<<endl;
    while(c8<=min)
      c8 = cm+sqrt(cs)*NR::gasdev(seed);
    d8 = exp(c8);
  }
  return d8;
}

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
  return 0.4*log(10.0)*phistar*pow(10.,-0.4*(M-Mstar)*(alpha+1))
    *exp(-1*pow(10.,-0.4*(M-Mstar)));
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
  for(double mag=-23.5; mag <=-16; mag+=0.0001){
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
  for(double mag=-23.5; mag <=-16; mag+=0.0001){
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
  ifstream file(filename.c_str());
  if (file.fail()) {
    cerr<<"error: cannot open "<<filename<<endl;
    exit(1);
  }
  else cout<<"reading "<<filename<<endl;
  int entries=3801;
  magnitude.resize(entries);
  lumnumdens.resize(entries);
  for (int ii=0; ii<entries; ii++) {
    file>>magnitude[ii]>>lumnumdens[ii];
  }
  file.close();
}


double LumNumberDensityInterp(double M){
  if((M<-25)||(M>-5)){
    cerr<<"Magnitude out of bounds in LumNumberDensity function: "<<M<<endl;
    exit(1);
  }
  vector<double>::iterator pos;
  pos=upper_bound(magnitude.begin(), magnitude.end(), M);
  int i = distance(magnitude.begin(),pos);
  double delta_x = magnitude[i] - magnitude[i-1];
  double dx      = M        - magnitude[i-1];
  double delta_y = lumnumdens[i] - lumnumdens[i-1];
  double y = dx/delta_x*delta_y+lumnumdens[i-1];

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

double NdensMagnitude(double ndens){
  double mmin = -23.;
  double mmax = -5;
  double diff = 1.;
  double nndensmid;
  double mmid;
  //cout<<ndens<<" "<<LumNumberDensity(mmin)<<" "<<LumNumberDensity(mmax)<<" ";
  while(diff>0.0001){
    mmid=(mmin+mmax)/2.;
    //    cout<<"range"<<ndens<<" "<<mmin<<" "<<mmax<<endl;
    double nndens1 = LumNumberDensityInterp(mmin);
    double nndens2 = LumNumberDensityInterp(mmax);
    nndensmid = LumNumberDensityInterp((mmin+mmax)/2);
    if((ndens<nndens1)||(ndens>nndens2)){
      cout<<mmin<<" "<<mmid<<" "<<mmax<<endl;
      cout<<"bad val"<<nndens1<<" "<<ndens<<" "<<nndens2<<endl;
    }
    assert((ndens>nndens1)&&(ndens<nndens2));
    diff = mmax-mmin;
    if(diff<0.0001) break;
    if(ndens>nndensmid) mmin=mmid;
    else mmax=mmid;
  }
  return mmid;
}


void ReadPDFs(vector <FiveTuple> &denspdf, vector <float> &mags){
  string filename = "denspdf/"+simlabel+"_denspdf.dat";
  ifstream file(filename.c_str());
  if (file.fail()) {
    cerr<<"error: cannot open "<<filename<<endl;
    exit(1);
  }
  else
    cerr<<"reading "<<filename<<endl;
 
  while(file){
    //float tmp;
    float mag, cm, fm, fs, p;
    //file>>mag>>tmp>>tmp>>cm>>fm>>fs>>p>>tmp>>tmp;
    file>>mag>>cm>>fm>>fs>>p;
    FiveTuple fTup(cm,fm,fs,p);
    mags.push_back(mag);
    denspdf.push_back(fTup);
  }

  for(int i=1;i<denspdf.size()-1;i++){
    float cm, fm, fs, p;
    //    denspdf[i].Print();
    cm = denspdf[i][0];
    fm = (denspdf[i-1][2]+denspdf[i][2]+denspdf[i+1][2])/3.;
    fs = (denspdf[i-1][3]+denspdf[i][3]+denspdf[i+1][3])/3.;
    p = (denspdf[i-1][4]+denspdf[i][4]+denspdf[i+1][4])/3.;
    FiveTuple fTup(cm,fm,fs,p); 
    denspdf[i] = fTup;
    //    denspdf[i].Print();
    //    for(int j=2;j<=4;j++){
      //      cout<<j<<" "<<denspdf[i-1][j]<<" "
      //  <<denspdf[i][j]<<" "
      //  <<denspdf[i+1][j]<<endl;
      //float newval = denspdf[i][j];
      
      //}
  }


  cout<<"Read"<<denspdf.size()<<"lines from denspdf";
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

  vector <FiveTuple> denspdf;
  vector <float> dmagbins;
  ReadPDFs(denspdf, dmagbins); 
  assert(denspdf.size()== dmagbins.size());
  float pdfbinsize = dmagbins[0]-dmagbins[1];
  int ind = dmagbins.size()-2;
  FiveTuple fTup = denspdf[ind];
  float cur_mag = dmagbins[ind];
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

  for(double mag=-23.7; mag <=-16; mag+=0.001){
    if(mag>cur_mag+0.5*pdfbinsize){
      if(ind >= 0){
	--ind;
	fTup = denspdf[ind];
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
      for(int i=0;i<diff;i++){
	Galaxy * galaxy = new Galaxy(mag,ngal,fTup.LocalDens());
 	galaxies.push_back(galaxy);
	ngal++;
      }
      if (galaxies.size()>=n) break;
  }  
  cout<<"[GetGalaxies] Done assigning galaxies "<<galaxies.size()<<endl;
  return galaxies;
}


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
                                                        
  for(double mag=-23.7; mag <=-16; mag+=0.001){

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

float SelectGalaxyZ()
{
  //selects a random radius for galaxy s.t. prob /propto r^3
  float rn = drand48();
  float rm = cosmo.RofZ(ZREDMAX);
  //rn *= rm*rm*rm*rm;
  //rn = pow(rn, 0.25);
  rn *= rm*rm*rm;
  rn = pow(rn, 1.0/3.0);
  float z = cosmo.ZofR(rn);
  return z;
}

float SelectGalaxyZ(float rmax)
{
  //selects a random radius for galaxy s.t. prob /propto r^3
  float rn = drand48();
  float rm = rmax;
  //rn *= rm*rm*rm*rm;
  //rn = pow(rn, 0.25);
  rn *= rm*rm*rm;
  rn = pow(rn, 1.0/3.0);
  float z = cosmo.ZofR(rn);
  return z;
}

vector <Galaxy *> GetGalaxies(double vol, ChainEl chel){
  integrator LF_integrator(&LF, 0, 1.0E-05, 0.01);
  cout<<vol<<endl;
  int n = (int) (LumNumberDensity(Magmin_pdf)*vol);//the 1.2 is padding
#ifdef MAG_LIMITED
  double dz = 0.01;
  double zNow = ZREDMIN;
  while(zNow < ZREDMAX){
    double r1 = cosmo.RofZ(zNow);
    double r2 = cosmo.RofZ(zNow+dz);
    double vol = 4.0/3.0*PI*(r2*r2*r2-r1*r1*r1);
    double dl = (r1+0.0001)*1e6*(1+cosmo.ZofR(r1))*(1+cosmo.ZofR(r1));
    double magmin_z = oMagMin - 5*(log10(dl)-1);
    cout<<" z = "<<zNow<<", magmin = "<<magmin_z<<", vol = "<<vol;
    int nnew = 0;
    if(magmin_z > Magmin)
      {
	cout<<", adding galaxies....";
	//nnew = (int) LF_integrator.Integrate(Magmin, magmin_z)*vol;
	nnew = (int) (LumNumberDensity(Magmin_pdf, magmin_z)*vol);
	n += nnew;
      }
    cout<<", new galaxies = "<<nnew<<endl;
    zNow += dz;
  }
#endif

  PRNT("GetGalaxies:", n);
  vector <Galaxy*> galaxies;
  galaxies.reserve(n);
  double volume = vol;
  cout<<volume<<" "<<sim.CubeVolume()<<endl;
  int ngal = 0;
  double LowestMag;

/*
#ifdef MAG_LIMITED
  //#ifdef DIM_GALAXIES
  for(double mag=-23.7; mag <0; mag+=0.001){
#else
  for(double mag=-23.7; mag <=-16; mag+=0.001){
#endif
      double num = LF_integrator.Integrate(-30, mag)*volume;
      int diff = ((int) floor(num))-galaxies.size();
      if (num>n) diff = n-galaxies.size();
      for(int i=0;i<diff;i++){
	float zGal = SelectGalaxyZ();
	float ThisMag = mag + Q*(zGal - 0.1);
#ifdef MAG_LIMITED
	float dl = cosmo.RofZ(zGal)*1e6*(1+zGal)*(1+zGal);
	float omag = ThisMag + 5*(log10(dl)-1);
	if (omag > oMagMin && mag > Magmin)
	  continue;
#endif
	float ThisMStar = Mstar + Q*(zGal - 0.1);
	float lum = pow(10.,-0.4*(ThisMag-ThisMStar));
	FiveTuple fTup(chel.cmean(),chel.fmean(lum),chel.fsig(lum),chel.ffrac(lum)); 
	Galaxy * galaxy = new Galaxy(ThisMag,ngal,fTup.LocalDens());
	galaxy->zGal(zGal);
 	galaxies.push_back(galaxy);
	ngal++;
      }
      LowestMag = mag;
      if (galaxies.size()>=n) break;
  }  
  */

  
  //First we do the main, volume-limited sample
  for(double mag=-23.7; mag <=Magmin; mag+=0.001){
    double num = LF_integrator.Integrate(-30, mag)*volume;
    int diff = ((int) floor(num))-galaxies.size();
    if (num>n) diff = n-galaxies.size();
    for(int i=0;i<diff;i++){
      float zGal = SelectGalaxyZ();
      float ThisMag = mag + Q*(zGal - 0.1);
      float ThisMStar = Mstar + Q*(zGal - 0.1);
      float lum = pow(10.,-0.4*(ThisMag-ThisMStar));
      FiveTuple fTup(chel.cmean(),chel.fmean(lum),chel.fsig(lum),chel.ffrac(lum)); 
      Galaxy * galaxy = new Galaxy(ThisMag,ngal,fTup.LocalDens());
      galaxy->zGal(zGal);
      galaxies.push_back(galaxy);
      ngal++;
    }
    LowestMag = mag;
    if (galaxies.size()>=n) break;
  }  
  
#ifdef MAG_LIMITED
  //Second the Magnitude-limited part
  cout<<"Volume Limited sample has "<<galaxies.size()<<" galaxies."<<endl;
  cout<<"Generating the Magnitude Limited Sample"<<endl;
  for(double mag=Magmin; mag <0; mag+=0.001){
    float DlMax = pow(10.0, 0.2*(oMagmin-mag)+1);
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
	float ThisMag = mag + Q*(zGal - 0.1);
	float dl = cosmo.RofZ(zGal)*1e6*(1+zGal)*(1+zGal);
	float omag = ThisMag + 5*(log10(dl)-1);
	if (omag > oMagMin && mag > Magmin)
	  continue;
	float ThisMStar = Mstar + Q*(zGal - 0.1);
	float lum = pow(10.,-0.4*(ThisMag-ThisMStar));
	FiveTuple fTup(chel.cmean(),chel.fmean(lum),chel.fsig(lum),chel.ffrac(lum)); 
	Galaxy * galaxy = new Galaxy(ThisMag,ngal,fTup.LocalDens());
	galaxy->zGal(zGal);
 	galaxies.push_back(galaxy);
	ngal++;
      }
      LowestMag = mag;
      if (galaxies.size()>=n) break;
  }  
#endif
  

  //  cout<<"[GetGalaxies] Done assigning galaxies "<<galaxies.size()<<endl;
  cout<<"[GetGalaxies] Done generating "<<galaxies.size()<<" down to magnitude "<<LowestMag<<endl;
  return galaxies;
}

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

