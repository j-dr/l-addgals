#include <vector>
#include <cassert> 
#include <unistd.h>     //for sleep()
#include <fstream>
#include <stdio.h>
#include <stdlib.h> //for qsort
#include <algorithm>
#include "galaxy.h"
#include "particle.h"
#include "halo.h"
#include "global_vars.h"
#include "ParameterDatabase.h"
#include "StringDatabase.h"
#include "ReadParameters.h"

extern double normal_random(float mean, float stddev);
extern double ranf(void);

/*
struct MrID{
  float Mr;
  int gid;
  float Z;
};
*/

class MrID{
public: 
  MrID(float Mr, int gid, float Z):magnitude(Mr), id(gid), Zred(Z){};
  float Mr()const{return magnitude;};
  int gid()const{return id;};
  float Z()const{return Zred;};
private:
  float magnitude;
  int id;
  float Zred;
};



//int CompareByMr (const struct MrID & a, const struct MrID & b)
int CompareByMr (MrID * a, MrID * b)
{
  return a->Mr() > b->Mr();
}

//int CompareByZ (const struct MrID & a, const struct MrID & b)
int CompareByZ (MrID * a, MrID * b)
{
  return a->Z() < b->Z();
}

// compare galaxy pointers by local mass density
bool DLessGal2(Galaxy * a, Galaxy * b)
{
  return a->Dist8() < b->Dist8(); 
}


#ifdef DENSPDF_IN_MAG_BIN
void recalculate_d8(vector <Galaxy *> &galaxies);
#endif

void read_growthfcn(vector <float> zarr, vector <float> darr){
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

//void SortGals(int &LLBins, int * &NInBin, int * &BinStart, struct MrID * GalMrID, int size, float GalZTol)
void SortGals(int &LLBins, int * &NInBin, int * &BinStart, vector <MrID *> &GalMrID, int size, float GalZTol)
{
  //cout<<"Before sort by z: "<<GalMrID[0]->Z()<<" "<<GalMrID[1]->Z()<<" "<<GalMrID[2]->Z()<<" ... "<<GalMrID[size-2]->Z()<<" "<<GalMrID[size-2]->Z()<<endl;
  sort(&GalMrID[0], &GalMrID[size], CompareByZ);
  //cout<<"Did initial sort by z: "<<GalMrID[0]->Z()<<" "<<GalMrID[1]->Z()<<" "<<GalMrID[2]->Z()<<" ... "<<GalMrID[size-2]->Z()<<" "<<GalMrID[size-2]->Z()<<endl;

  //Determine the numbeer of bins
  LLBins = 0;
  float tz = 0.0;
  while (tz < ZREDMAX){
    tz += (1+tz)*GalZTol;
    LLBins++;
  }
  NInBin = new int[LLBins];
  BinStart = new int[LLBins];
  for(int i=0;i<LLBins;i++){
    NInBin[i] = 0;
    BinStart[i] = 0;
  }

  //cout<<" going to sort "<<LLBins<<" bins in order of magnitude"<<endl;
  //Sort the particles in the z-bins in order of magnitude
  int FirstGal, LastGal;
  float z0, z1;
  FirstGal = 0;
  z0 = 0.0;
  for(int ib=0;ib<LLBins;ib++){
    BinStart[ib] = FirstGal;
    z1 = z0+(1+z0)*GalZTol;
    int ig = FirstGal;
    while(ig < size && GalMrID[ig]->Z() < z1){
      ig++;
      if(ig >= size)
	break;
    }
    LastGal = ig;
    NInBin[ib] = LastGal - FirstGal;
    //cout<<" going to sort bin "<<ib<<" with first gal = "<<FirstGal<<", last gal = "<<LastGal<<" out of a total of "<<size<<" galaxies (galaxies.size = "<<GalMrID.size()<<")"<<endl;
    if (NInBin[ib] > 0)
      {
	//cout<<" sorting bin..."<<endl;
	sort(&GalMrID[FirstGal], &GalMrID[LastGal], CompareByMr);
	//cout<<"               done"<<endl;
      }
    //cout<<"The sorted Mr distribution: "<<GalMrID[0]->Mr()<<" "<<GalMrID[1]->Mr()<<" "<<GalMrID[2]->Mr()<<" ... "<<GalMrID[size-2]->Mr()<<" "<<GalMrID[size-1]->Mr()<<endl;
    //cout<<"The (should be screwed up) z-distribution: "<<GalMrID[0]->Z()<<" "<<GalMrID[1]->Z()<<" "<<GalMrID[2]->Z()<<" ... "<<GalMrID[size-2]->Z()<<" "<<GalMrID[size-2]->Z()<<endl;
    z0 = z1;
    FirstGal = LastGal;
  }

}

//int SelectGalaxy(int &LLBins, int * &NInBin, int * &BinStart, struct MrID * GalMrID, float mr, float zHalo, float zGalTol, vector <Galaxy *> &galaxies)
int SelectGalaxy(int &LLBins, int * &NInBin, int * &BinStart, vector <MrID *> &GalMrID, float mr, float zHalo, float zGalTol, vector <Galaxy *> &galaxies)
{
  int Bin = 0;
  float z = 0.0;
  while (z < zHalo){
    Bin++;
    z += (1+z)*zGalTol;
  }
  Bin--;
  if (Bin >= LLBins){
    cout<<"Error!  zBin is too big!  zBin = "<<Bin<<" z = "<<zHalo<<" LLBins = "<<LLBins<<endl;
  }
  //First we do a coarse search to get the general location
  int iGalSelect;
  for(iGalSelect=BinStart[Bin];iGalSelect<BinStart[Bin]+NInBin[Bin];iGalSelect+= 1000){
    //check direction of inequality -- search until we find something brighter
    if(GalMrID[iGalSelect]->Mr() < mr)
      break;
  }
  //Now we actually find the best-fit halo
  int iGalStart = iGalSelect - 1000;
  if (iGalStart < BinStart[Bin])
    iGalStart = BinStart[Bin];
  for (iGalSelect=iGalStart;iGalSelect<BinStart[Bin]+NInBin[Bin];iGalSelect++){
    if(GalMrID[iGalSelect]->Mr() <= mr && !(galaxies[GalMrID[iGalSelect]->gid()]->Central()))
      break;
  }
  if(iGalSelect == BinStart[Bin]+NInBin[Bin])
    iGalSelect = -1;

  return iGalSelect;
}
      
  


//We add the BCGs by looping through our halos and, if the halo 
//is in the volume, we select a galaxies from our LF too assign
//if we can find a good match in mr, else we create a new galaxy.
//If the halo is not in our volume, we just create a new galaxy.

void add_bcgs(vector <Particle *> &particles, vector <Galaxy *> &galaxies, vector  <Halo *> &halos)
{
  //This sort matches up with a sort done in assigngals
  //sort(galaxies.begin(),galaxies.end(),DLessGal2);
  //cout<<"Putting in the BCGS"<<endl;

  string outcheck="./BCG_check.dat";
  ofstream outcheckfile(outcheck.c_str());

  //string galcheck="./Galaxy_check.dat";
  //ofstream galcheckfile(galcheck.c_str());
  //for(int ig=0;ig<galaxies.size();ig++){
  //  galcheckfile<<galaxies[ig]->Mr()<<endl;
  //}

  ///read the growth fcn so we can resize the galaxies
  vector <float> zarr;
  vector <float> darr;
  read_growthfcn(zarr, darr);

  int ExtraGals = 0;    //The number of new galaxies we've added
  int ExtraGalsInVol = 0;    //The number of new galaxies we've added to our volume

  //Do a LL sort to get galay luminosity sorted but not change the structure
  //MrID GalMrID[galaxies.size()];
  cout<<" Doing Sort for BCGs"<<endl;
  cout<<"   (ngals = "<<galaxies.size()<<")"<<endl;
  //  struct MrID GalMrID[galaxies.size()];
  vector <MrID *> GalMrID;
  GalMrID.reserve(galaxies.size());
  cout<<"   (Made the structure for id's)"<<endl;
  cout<<" Initializing Mr's"<<endl;
  for(int ig=0;ig<galaxies.size();ig++){
    //GalMrID[ig].Mr = galaxies[ig]->Mr() + Q*(galaxies[ig]->zGal()-0.1);
    //    GalMrID[ig].Mr = galaxies[ig]->Mr();
    //    GalMrID[ig].gid = ig;
    //    GalMrID[ig].Z = galaxies[ig]->zGal();
    //GalMrID[ig]->Mr = galaxies[ig]->Mr();
    MrID * tGalMrID = new MrID(galaxies[ig]->Mr(), ig, galaxies[ig]->zGal());
    GalMrID.push_back(tGalMrID);
  }
  cout<<" Initialized arrays"<<endl;
  int Nsort = galaxies.size();
  int ElementSize = sizeof(struct MrID);
  cout<<"  Starting Sort."<<endl;
  //sort(&GalMrID[0], &GalMrID[galaxies.size()-1], CompareByMr);
  //SortGals(&GalMrID, galaxies.size()-1);

  //Why did I earlier set GalZTol separately? ->I think to prevent small number statistics from causing us to add too many bright galaxies
  float GalZTol = zTol;
  if (GalZTol < 0.01)
    GalZTol = 0.01;

  int LLBins = 0;
  int * NInBin;
  int * BinStart;
  SortGals(LLBins, NInBin, BinStart, GalMrID, galaxies.size()-1, GalZTol);
  //sort(&GalMrID[0], &GalMrID[galaxies.size()-1], CompareByZ);
  cout<<"finished sort"<<endl;
  //cout<<GalMrID[0]->Mr()<<" "<<GalMrID[1]->Mr()<<" "<<GalMrID[2]->Mr()<<" ... "<<GalMrID[galaxies.size()-2]->Mr()<<" "<<GalMrID[galaxies.size()-1]->Mr()<<endl;

  float NextPercent = 0.0;
  int OriginalGalSize = galaxies.size();   // # of galaxies before we need to add any BCGs
  cout<<" Inserting BCGs into "<<halos.size()<<" halos."<<endl;
  for(int hi=0; hi<halos.size();hi++){
    if (float(hi)/float(halos.size()-1) >= NextPercent){
      int np = int(NextPercent*100);
      cout<<"  "<<np<<"% done"<<endl;
      NextPercent += 0.1;
    }
    int pid = 0;
    
    if(!(halos[hi]->InVol()))
      continue;

    float this_growth_fcn = 1.0;
    for (int i=0;i<zarr.size();i++){
      this_growth_fcn = darr[i];
      if (zarr[i] >= halos[hi]->ZredReal())
	continue;
    }
	       
    double m200 = halos[hi]->M();
    if (m200 < BCG_Mass_lim)
      continue;
    //Simple Power Law
    //double m14 = m200/1e14;
    //double lum = 4*pow(m14,0.3);
    //Zeng, Coil, Zehavi
    float Mc = 3.7e9;
    float a = 29.78;
    float b = 29.5;
    float k = 0.0255;
    float L0 = 2.8e9;
    L0 *= (5./8.);  // correction fudge factor to get Sarah's plot correct
    //float L0 = 3.36e9; //recalculated to match sham
    L0 *= 7.80e-11;  //Convert from solar to L*
    double loglum = log10(L0) + a*log10(m200/Mc) - (1./k)*log10(1.+pow(m200/Mc,b*k));
    double lum = pow(10., loglum);
    double lum_before = lum;
    lum = pow(10.0,normal_random(log10(lum),0.15));
    double mr = -2.5*log10(lum) + Mstar;


    /*
    //now we fortget about the above and do our simple power-law
    double fit_00 = -10.8998;
    double fit_01 = 1.17409;
    double fit_10 = -0.818732;
    double fit_11 =  -0.115187;

    double fit0 = fit_00 + fit_01*halos[hi]->Zred();
    double fit1 = fit_10 + fit_11*halos[hi]->Zred();

#ifdef SNAPSHOT
    fit0 = fit_00;
    fit1 = fit_10;
#endif

    mr = fit0 + fit1*log10(m200);
    lum = pow(10.0, -0.4*(mr - Mstar));
    lum_before = lum;
    lum = pow(10.0, normal_random(log10(lum), 0.15));
    mr = -2.5*log10(lum) + Mstar;
    */

    //We de-evolve the galaxy to keep the catalog consistent
    /*
#ifdef SNAPSHOT
    if (mr > Magmin) 
#else
    if (evolve_mag(mr, halos[hi]->Zred()) > Magmin) 
#endif
    */    
    if (mr > Magmin) 
      continue;

    //find a galaxy with this magnitude
    int iGalSelect;
    int NewFlag = 0;
    //float zCut = ztol*(1+halos[hi]->Zred())
    
    if(halos[hi]->InVol()){
      iGalSelect = SelectGalaxy(LLBins, NInBin, BinStart, GalMrID, mr, halos[hi]->Zred(), GalZTol, galaxies);

      //if (iGalSelect == -1)
      //cout<<"iGalSelect = -1 for halo "<<hi<<"!  Desired Mr = "<<mr<<endl;


      //if not a galaxy bright enought, we just make one
      if(iGalSelect == -1){
	Galaxy * galaxy = new Galaxy(mr,galaxies.size(),1e12);
	galaxies.push_back(galaxy);
	ExtraGals++;
	ExtraGalsInVol++;
	iGalSelect = galaxies.size()-1;
	NewFlag = 1;
      }


    }
    else{//The halo is outside the volume, we we may have a particle linking to it anyway so create a new galaxy
      Galaxy * galaxy = new Galaxy(mr,galaxies.size(),1e12);
      galaxies.push_back(galaxy);
      ExtraGals++;
      iGalSelect = galaxies.size()-1;
      NewFlag = 1;
    }
      
    //create a new particle for the halo
    float xfac = 1./(sim.LengthUnit());
    Point xx(halos[hi]->X()*xfac,halos[hi]->Y()*xfac,halos[hi]->Z()*xfac);
    Point vv(halos[hi]->Vx(),halos[hi]->Vy(),halos[hi]->Vz());
    float dist8 = 1e12;  //Just a hack -- Low density to put them at the end of the list
    Particle * particle = new Particle(xx,vv,dist8);
    particle->Hid(hi);
    particles.push_back(particle);
    if(NewFlag){
      galaxies[iGalSelect]->P(particles[particles.size()-1]);
      galaxies[iGalSelect]->DefineCentral();
      if(halos[hi]->InVol())
	outcheckfile<<m200<<" "<<lum_before<<" "<<lum<<" "<<mr<<" "<<galaxies[iGalSelect]->Mr()<<" "<<halos[hi]->Zred()<<" 1 "<<endl;
      halos[hi]->Brightest(iGalSelect);
      halos[hi]->Central(iGalSelect);
      halos[hi]->CentralDistance(0.);
      particles[particles.size()-1]->MakeGal(iGalSelect);
    } else {
      galaxies[GalMrID[iGalSelect]->gid()]->P(particles[particles.size()-1]);
      galaxies[GalMrID[iGalSelect]->gid()]->DefineCentral();
      if(halos[hi]->InVol())
	outcheckfile<<m200<<" "<<lum_before<<" "<<lum<<" "<<mr<<" "<<galaxies[GalMrID[iGalSelect]->gid()]->Mr()<<" "<<halos[hi]->Zred()<<" 0 "<<endl;
      halos[hi]->Brightest(GalMrID[iGalSelect]->gid());
      halos[hi]->Central(GalMrID[iGalSelect]->gid());
      halos[hi]->CentralDistance(0.);
      particles[particles.size()-1]->MakeGal(GalMrID[iGalSelect]->gid());
    }
    
  }

  cout<<"Had to add "<<ExtraGals<<" extra galaxies for BCGs"<<endl;
  cout<<"  "<<ExtraGalsInVol<<" of which where in the specified volume"<<endl;
  int NCentralsFound = 0;
  for(int igc = 0;igc<galaxies.size();igc++){
    if (galaxies[igc]->Central())
      {
	int hi = galaxies[igc]->P()->Hid();
	double halor = galaxies[igc]->P()->Distance(halos[hi]->Position());
	//cout<<"a BCG halor = "<<halor<<endl;
	NCentralsFound++;
      }
  }
  cout<<"Confirmed the existance of "<<NCentralsFound<<" centrals."<<endl;

#ifdef DENSPDF_IN_MAG_BIN
  cout<<"Recalculating d8 for the non-BCGs..."<<endl;
  recalculate_d8(galaxies);
  cout<<"                                     done"<<endl;
#endif
}
