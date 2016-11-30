#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <list>
#include <strstream>
#include <iostream>
#include <sstream>
#include <string>

#include "galaxy.h"
#include "biniostream.h"
#include "point.h"
#include "particle.h"
#include "halo.h"
#include "singleton.h"
#include "functions.h"
#include "constants.h"
//#include "galaxy.h"
#include "box.h"

vector <Galaxy *> ReadSHAMGalaxies(string shamfile)
{
  //ifstream rfile(rockfile.c_str());
  ifstream sfile(shamfile.c_str());

  //if (rfile.fail()){
  //  cerr<<"Error opening rockstar file: "<<rockfile<<endl;
  //  exit(1);
  //}
  if (sfile.fail()){
    cerr<<"Error opening sham file: "<<shamfile<<endl;
    exit(1);
  }

  float mr, d8;
  int id;
  string tmps;
 
  sfile>>tmps;
  sfile>>tmps;
  sfile>>tmps;
  int ng = 0;
  vector <Galaxy*> galaxies;
  sfile>>id>>mr>>d8;
  while(!sfile.eof()){
    //sfile>>id>>mr>>d8;
    if (ng < 10) cout<<id<<" "<<mr<<" "<<d8<<endl;
    if (mr <= Magmin){
      Galaxy *galaxy = new Galaxy(mr, id, d8);
      galaxies.push_back(galaxy);
      galaxies[ng]->zGal(0.0);
      ng++;
    }
    sfile>>id>>mr>>d8;
  }
  sfile.close();

  //output some diagnostics of our galaxy distribution
  cout<<"[ReadSHAMGalaxies] Done getting galaxy magnitudes "<<galaxies.size()<<endl;
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

void get_mhost(vector <Halo *> &halos, vector <Galaxy *> &galaxies){

  cout<<"starting get_mhost"<<endl;
  //make our lookup table to get HOD information
  int max_hid = 0;
  for (int ih=0;ih<halos.size();ih++)
    if (halos[ih]->Id() > max_hid) max_hid = halos[ih]->Id();
  int* lookuptable = new int[max_hid+1];
  for(int ih=0;ih<halos.size();ih++)
    lookuptable[halos[ih]->Id()] = ih;

  //loop through galaxies to find host halo
  for(int ih=0;ih<halos.size();ih++){
    int pid = lookuptable[halos[ih]->Host()];
    while (lookuptable[halos[pid]->Host()] != pid)
      pid = lookuptable[halos[pid]->Host()];
    galaxies[ih]->Mhost(halos[pid]->M());
    //if (ih < 10) cout<<"Mhost value: "<<galaxies[ih]->Mhost()<<endl;
    galaxies[ih]->H(halos[ih]);
  }
}

#ifdef JUST_COLORS

void read_galaxies(vector <Particle *> &particles, vector <Galaxy *> &galaxies, vector <Halo *> &halos){
  string filename;

  filename = galaxy_file;

  ifstream file(filename.c_str());
  if (file.fail()){
    cerr<<"Error reading galaxies:  Can't open file "<<filename<<endl;
    exit(1);
  } else {
    cout<<"Reading galaxies from file "<<filename<<endl;
  }

  float px, py, pz, vx, vy, vz, mr, d8, m200, r200;
  int hid, gid, cent;
  string tmps;
#ifdef BOLSHOI
  float vmax, vpeak, mpeak, zpeak, mag0, mag1, mag2, mag3, mag4, mag5, mag6;
  string tmps;
  float tgid, thid;
  while(tmps != "0.4")
    {
      file>>tmps;
      cout<<tmps<<endl;
    }
  //file>>tmps;
  //cout<<tmps<<endl;
#elif MASSIVE_BLACK
  file>>tmps>>tmps>>tmps>>tmps>>tmps;
  cout<<tmps<<endl;
#endif

  int ngal = 0;
  int nhalo = 0;
  float zmin = 100.;
  float maxz = 0.;
  float ramin = 500.;
  float ramax = 0.;
  float decmin = 500.;
  float decmax = 0.;
  while(file){
#ifdef BOLSHOI
    file>>tgid>>thid>>px>>py>>pz>>vx>>vy>>vz>>m200>>r200>>vmax>>vpeak>>mpeak>>zpeak>>mag0>>mag1>>mag2>>mag3>>mag4>>mag5>>mag6;
    hid = int(thid);
    gid = int(tgid);
    //out<<gid<<" "<<px<<endl;
    mr = mag3;
    d8 = 0;
    cent = 0;
    if (hid == 0)
      {
	cent = 1;
	hid = gid;
      }
#elif MASSIVE_BLACK
    file>>gid>>px>>py>>pz>>mr;
    hid = int(gid);
    d8 = 0;
    cent = 0;
    vx = 0.0;
    vy = 0.0;
    vz = 0.0;
    m200 = 0.0;
#else
    file>>px>>py>>pz>>vx>>vy>>vz>>mr>>hid>>m200>>r200>>gid>>d8>>cent;
#endif
    Point xx(px/sim.LengthUnit(),py/sim.LengthUnit(),pz/sim.LengthUnit());
    Point vv(vx,vy,vz);

    Particle * particle = new Particle(xx,vv,d8);
    particle->SetHaloId(-1);
    //Galaxy * galaxy = new Galaxy(mr,ngal,d8);
    Galaxy * galaxy = new Galaxy(mr,gid,d8);

    if (cent && m200 > BCG_Mass_lim) {
      galaxy->DefineCentral();
      int np = 0;
      Halo * halo = new Halo(xx,vv,m200,particle->Zred(),np,hid,0.,r200);
      halos.push_back(halo);
      particle->SetHaloId(nhalo);
      nhalo++;
    }
    galaxies.push_back(galaxy);
    particles.push_back(particle);
    galaxies[ngal]->zGal(particles[ngal]->Zred());
    galaxies[ngal]->P(particles[ngal]);
    particles[ngal]->MakeGal(ngal);
    //we need to de-evolve the galaxy magnitude
    //DeEvolveGal(galaxies[ngal]);

#ifdef SNAPSHOT
    if(particle->Z()>zmax->GetVal())
      zmax->SetVal(particle->Z());
#else
    if(particle->Zred()>zmax->GetVal())
      zmax->SetVal(particle->Zred());
#endif
    ///*
    if (ngal < 10){
      particles[ngal]->PosPrint();
      cout<<particles[ngal]->Ra()<<" "<<particles[ngal]->Dec()<<" "<<particles[ngal]->Zred()<<endl;
      cout<<galaxies[ngal]->Ra()<<" "<<galaxies[ngal]->Dec()<<" "<<galaxies[ngal]->Z()<<" "<<galaxies[ngal]->Zbin()<<endl;
      //particles[ngal]->Write();
    }
    //*/
    ngal++;
  }

  //delete the last element which is just the last line read in twice
  galaxies.erase(galaxies.end()-2, galaxies.end()-1);
  particles.erase(particles.end()-2, particles.end()-1);

  //Just re-check the range of particles/galaxies
  for(int i=0;i<galaxies.size();i++){
    if (galaxies[i]->Z() > maxz)
      maxz = galaxies[i]->Z();
    if (galaxies[i]->Z() < zmin)
      zmin = galaxies[i]->Z();
    if (galaxies[i]->Ra() > ramax)
      ramax = galaxies[i]->Ra();
    if (galaxies[i]->Ra() < ramin)
      ramin = galaxies[i]->Ra();
    if (galaxies[i]->Dec() > decmax)
      decmax = galaxies[i]->Dec();
    if (galaxies[i]->Dec() < decmin)
      decmin = galaxies[i]->Dec();
  }

  cout<<"Finished reading "<<galaxies.size()<<" galaxies."<<endl;
  cout<<"  min/max z = "<<zmin<<" / "<<maxz<<endl;
  cout<<"  min/max ra = "<<ramin<<" / "<<ramax<<endl;
  cout<<"  min/max dec = "<<decmin<<" / "<<decmax<<endl;
  cout<<"First 3 galaxy positions:"<<endl;
  cout<<"  "<<" "<<galaxies[0]->P()->X()<<" "<<galaxies[0]->P()->Y()<<" "<<galaxies[0]->P()->Z()<<endl;
  cout<<"  "<<" "<<galaxies[1]->P()->X()<<" "<<galaxies[1]->P()->Y()<<" "<<galaxies[1]->P()->Z()<<endl;
  cout<<"  "<<" "<<galaxies[2]->P()->X()<<" "<<galaxies[2]->P()->Y()<<" "<<galaxies[2]->P()->Z()<<endl;
  cout<<"Last 3 galaxy positions:"<<endl;
  int ig = galaxies.size()-3;
  cout<<"  "<<" "<<galaxies[ig]->P()->X()<<" "<<galaxies[ig]->P()->Y()<<" "<<galaxies[ig]->P()->Z()<<endl;
  ig = galaxies.size()-2;
  cout<<"  "<<" "<<galaxies[ig]->P()->X()<<" "<<galaxies[ig]->P()->Y()<<" "<<galaxies[ig]->P()->Z()<<endl;
  ig = galaxies.size()-1;
  cout<<"  "<<" "<<galaxies[ig]->P()->X()<<" "<<galaxies[ig]->P()->Y()<<" "<<galaxies[ig]->P()->Z()<<endl;

}

#endif
