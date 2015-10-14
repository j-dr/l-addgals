#include <vector>
#include <cassert> 
#include <unistd.h>     //for sleep()
#include <fstream>
//#include "hv2.h"
//#include "hv.h"
#include "galaxy.h"

int GalaxyZBin(float zRed);

void PrintMinMaxDens(vector <Galaxy *> &galaxies, vector <Particle *> &particles)
{
  float min_p_dens = particles[0]->Dist8();
  float max_p_dens = particles[0]->Dist8();
  float min_g_dens = galaxies[0]->Dist8();
  float max_g_dens = galaxies[0]->Dist8();
  int min_p = 0;
  int max_p = 0;
  int min_g = 0;
  int max_g = 0;
  for(int pi=1;pi<particles.size();pi++)
    {
      if (particles[pi]->Dist8() < min_p_dens)
	{
	  min_p_dens = particles[pi]->Dist8();
	  min_p = pi;
	}
      if (particles[pi]->Dist8() > max_p_dens)
	{
	  max_p_dens = particles[pi]->Dist8();
	  max_p = pi;
	}
    }
  for(int gi=1;gi<galaxies.size();gi++)
    {
      if (galaxies[gi]->Dist8() < min_g_dens)
	{
	  min_g_dens = galaxies[gi]->Dist8();
	  min_g = gi;
	}
      if (galaxies[gi]->Dist8() > max_g_dens)
	{
	  max_g_dens = galaxies[gi]->Dist8();
	  max_g = gi;
	}
    }
  cout<<"Min/Max particle densities = "<<min_p_dens<<", "<<max_p_dens<<endl;
  cout<<"For particles "<<min_p<<", "<<max_p<<endl;
  cout<<"Min/Max galaxies densities = "<<min_g_dens<<", "<<max_g_dens<<endl;
  cout<<"For galaxies "<<min_g<<", "<<max_g<<endl;
}


void print_galaxies(vector <Galaxy *> &galaxies, vector <Particle *> &particles, vector <Halo *> &halos, vector <GalSED> galseds, vector <int> sed_ids, vector <float> nndist, vector <float> nndist_percent, string outpfn, string outdfn, string outgfn, string outghfn, string outgzfn, string outrfn)
//void print_galaxies(vector <Galaxy *> &galaxies, vector <Particle *> &particles, vector <Halo *> &halos, vector <GalSED> galseds, vector <int> sed_ids, vector <float> nndist)
{
  ofstream outpfile(outpfn.c_str());
  ofstream outdfile(outdfn.c_str());
  ofstream outgfile(outgfn.c_str());
  ofstream outghfile(outghfn.c_str());
  ofstream outgzfile(outgzfn.c_str());
  ofstream outrfile(outrfn.c_str());

  int NNotPrinted = 0;
  int NNotPrintedInVol = 0;
  for(int gi=0;gi<galaxies.size();gi++){
    Galaxy * gal = galaxies[gi];
    Particle * p = gal->P();
    assert(p);  //this better be true since you removed the other ones.
    int hid = p->Hid();
#ifdef COLORS
    int sedid = sed_ids[gi];
    if((nndist[gi]>0)&&(sedid>=0)){
      if(galseds[sed_ids[gi]].CatId() == 773)
	cout<<"Missing galaxy info = "<<gal->Mr()<<endl;
      //      outgfile<<sed_ids[gi]<<" ";
      outgfile<<galseds[sed_ids[gi]].CatId()<<" ";
      //      SEDs[sed_ids[gi]].Write(outgfile);
#else //make a dummy for the color index
    if(1){
      outgfile<<0<<" ";
#endif
      gal->Write(outgfile);
      outgzfile<<gal->zGal()<<" "<<GalaxyZBin(gal->zGal())<<" "<<gal->Central()<<endl;
      //cout<<"Writing P info..."<<endl;
      p->Write(outpfile);
      //cout<<"Writing halo info..."<<endl;
      //cout<<" hid = ";
      //cout<<hid<<endl;
      if(hid<0 || hid>=halos.size()){
	//cout<<" hid = 0"<<endl;
	outghfile<<"0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "<<endl;
      }
      else{
	/*
	cout<<" Mass: "<<endl;
	cout<<halos[hid]->M()<<endl;
	cout<<" Ngal: ";
	cout<<halos[hid]->Ngal()<<endl;
	cout<<" r200: ";
	cout<<halos[hid]->R200()<<endl;
	cout<<" dist: ";
	cout<<p->Distance(halos[hid]->Position())<<endl;
	cout<<" Sig: ";
	cout<<halos[hid]->Sig()<<endl;
	cout<<" X: ";
	cout<<halos[hid]->X()<<endl;
	cout<<" Y: ";
	cout<<halos[hid]->Y()<<endl;
	cout<<" Z: ";
	cout<<halos[hid]->Z()<<endl;
	cout<<" Vx: ";
	cout<<halos[hid]->Vx()<<endl;
	cout<<" Vy: ";
	cout<<halos[hid]->Vy()<<endl;
	cout<<" Vz: ";
	cout<<halos[hid]->Vz()<<endl;
	cout<<" Dec: ";
	cout<<halos[hid]->Ra()<<endl;
	cout<<" Ra: ";
	cout<<halos[hid]->Dec()<<endl;
	cout<<" Zred: ";
	cout<<halos[hid]->Zred()<<endl;
	*/
	outghfile//<<in_vol<<" "
	  //<<hid<<" "
	  <<halos[hid]->Id()<<" "
	  <<halos[hid]->M()<<" "
	  <<halos[hid]->Ngal()<<" "
	  <<halos[hid]->R200()<<" "
	  <<p->Distance(halos[hid]->Position())<<" "
	  <<halos[hid]->Sig()<<" "
	  <<halos[hid]->X()<<" "
	  <<halos[hid]->Y()<<" "
	  <<halos[hid]->Z()<<" "
	  <<halos[hid]->Vx()<<" "
	  <<halos[hid]->Vy()<<" "
	  <<halos[hid]->Vz()<<" "
	  <<halos[hid]->Ra()<<" "
	  <<halos[hid]->Dec()<<" "
	  <<halos[hid]->Zred()<<endl;
	}
      //cout<<"Writing d8 info..."<<endl;
      //cout<<" gal d8: ";
      //cout<<gal->Dist8()<<endl;
      //cout<<" P d8: ";
      //cout<<gal->P()->Dist8()<<endl;
      outdfile<<gal->Dist8()<<" "<<gal->P()->Dist8()<<endl;
#ifdef COLORS
#ifdef COLORS_FROM_RELATIVE_DENSITY
      //cout<<"Writing nndist info..."<<endl;
      //cout<<" nndist: ";
      //cout<<nndist[gi]<<endl;
      //cout<<nndist_percent[gi]<<endl;
      outrfile<<nndist[gi]<<" "<<nndist_percent[gi]<<endl;
#else
      outrfile<<nndist[gi]<<endl;
#endif
#endif
    }
    else{
      NNotPrinted++;
      if(p->Save()){
	NNotPrintedInVol++;
	cout<<"[hv] didn't print galaxy "
	    <<gi<<" "
#ifdef COLORS
	    <<nndist[gi]<<" "
	    <<sedid<<" "
#endif
	    <<p->X()<<" "
	    <<p->Y()<<" "
	    <<p->Z()<<endl;
      }
    }
  }
    cout<<"Didn't print "<<NNotPrinted<<" galaxies.  Of these "<<NNotPrintedInVol<<" were in the specified volume."<<endl;  
}


