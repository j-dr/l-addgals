#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <list>
#include <strstream>
#include <iostream>
#include <sstream>
#include <string>

#include "biniostream.h"
#include "point.h"
#include "particle.h"
#include "halo.h"
#include "singleton.h"
#include "functions.h"


static const int MMAXINT=32767;


//
// reads all halos from file and returns vector of halo pointers
//
vector <Halo*> ReadMGSHalos(void){
  vector <Halo*> halos;

  string filename = datadir+"groups/grps.addgals.tot.dat";
  float mpart = 1.4224436e10;

  int hid = 0;

  ifstream file(filename.c_str());
  if (file.fail()) {
    cerr<<"[read_MGS error]: cannot open "<<filename<<endl;
    exit(1);
  } else {
    cout<<"[read_MGS]: opened file "<<filename<<" for reading halo information."<<endl;
  }
  while(file){
    float mtot, mdm, xMpc, yMpc, zMpc, vx, vy, vz, rdelta, junk1, junk2, junk3, siglos;
    int id, snap_num, ip;
    file>>id>>snap_num>>mtot>>mdm>>xMpc>>yMpc>>zMpc>>vx>>vy>>vz>>rdelta>>junk1>>junk2>>junk3;

    if (hid == 0){
      cout<<"First Halo Position = "<<xMpc<<" "<<yMpc<<" "<<zMpc<<endl;
    }
    
    siglos = 0.;
    Point vel(vx,vy,vz);   
    float xx, yy, zz;
    
    xx = xMpc;
    yy = yMpc;
    zz = zMpc;
    
    ip = int(mtot/mpart);
    float zred = cosmo.ZofR(sqrt(xMpc*xMpc + yMpc*yMpc + zMpc*zMpc));
    Point pos(xx, yy, zz);
    Halo * halo = new Halo(pos, vel, mtot, zred, ip, hid, siglos, rdelta);
    
    halos.push_back(halo);
    hid++;
      
  }
  file.close();
  cout<<"[read_MGS] Read "<<halos.size()<<" halos."<<endl;
  return halos;
}



int ReadMGS(int p, int t, vector <Particle *> &particles){
  int buf;
  unsigned int np = 0;
  float ramax = -1000.;
  float decmax = -1000.;
  float z_max = -1000.;
  float ramin = 1000.;
  float decmin = 1000.;
  float z_min = 1000.;

  int nbad = 352680;
  nbad = 1000000000;

  string phi1 = MakeString(p,2);
  string phi2 = MakeString(p+1,2);
  string theta1 = MakeString(t,2);
  string theta2 = MakeString(t+1,2);
  string FirstDir = datadir+"data/run.c1.s69.p"+phi1+"."+phi2+".t"+theta1+"."+theta2+".p0.d1/";
  string SecondDir = FirstDir+"run.s069.p"+phi1+".t"+theta1+".p0/";
  for(int ss=0;ss<160;ss++){
    string fstrp = SecondDir + "lc.pos.s"+MakeString(ss,3)+".p"+phi1+".t"+theta1+".p0.dat";
    string fstrv = SecondDir + "lc.vel.s"+MakeString(ss,3)+".p"+phi1+".t"+theta1+".p0.dat";
    string fstrh = SecondDir + "rnn.pos3.s"+MakeString(ss,3)+".p"+phi1+".t"+theta1+".p0.dat";
    string fstrhid = SecondDir + "hid.ids.s"+MakeString(ss,3)+".p"+phi1+".t"+theta1+".p0.dat";
  
    binifstream pfile(fstrp.c_str());
    binifstream vfile(fstrv.c_str());
    binifstream hfile(fstrh.c_str());
    ifstream hidfile(fstrhid.c_str());

    //confirm that files open
    if (pfile.fail()) {
      //cerr<<"error: cannot open file '" <<fstrp<<"'"<<endl;
      continue;
    }
    if (vfile.fail()) {
      cerr<<"error: cannot open file '" <<fstrv<<"'"<<endl;
      continue;
    }
    if (hfile.fail()) {
      cerr<<"error: cannot open file '" <<fstrh<<"'"<<endl;
      continue;
    }
#ifdef CATALOG
    if (hidfile.fail()) {
      //cerr<<"error: cannot open file '" <<fstrhid<<"'"<<endl;
      //continue;
    }
#endif

    int buf, np_file;
    pfile>>buf>>np_file>>buf>>buf;
    vfile>>buf>>np_file>>buf>>buf;
    hfile>>buf>>np_file>>buf>>buf;
    cout<<"reading pos file '"<<fstrp<<"' with "<<np_file<<" particles."<<endl;
#ifdef CATALOG
    //hidfile>>buf>>np_file>>buf>>buf;
#endif
    
    float xfac = 1./(sim.LengthUnit());
    for(int i=0;i<np_file;i++){
      float px,py,pz,vx,vy,vz,tmp;
      pfile>>px>>py>>pz;
      //zpos[i] = tmp;
      vfile>>vx>>vy>>vz;
      //zvel[i] = tmp;
      Point xx(px*xfac,py*xfac,pz*xfac);
      Point vv(vx,vy,vz);
      hfile>>tmp;
      Particle * particle = new Particle(xx,vv,tmp);
      //particle->PosPrint();
      if (particle->Ra() > ramax)
	ramax = particle->Ra();
      if (particle->Ra() < ramin)
	ramin = particle->Ra();
      if (particle->Dec() > decmax)
	decmax = particle->Dec();
      if (particle->Dec() < decmin)
	decmin = particle->Dec();
      if (particle->Zred() > z_max)
	z_max = particle->Zred();
      if (particle->Zred() < z_min)
	z_min = particle->Zred();

#ifdef CATALOG
      particle->Hid(0);
      if(particle->Save()){
	if (np > nbad)
	  cout<<" going to try and save "<<endl;
	particles.push_back(particle);
	np++;
	if(particle->Zred()>zmax->GetVal())
	  zmax->SetVal(particle->Zred());
	if (np > nbad)
	  cout<<"              successfully saved "<<np<<endl;
      }
      else{
	delete particle;
      }
#endif
#ifdef HALOASSIGN
      particles.push_back(particle);
      np++;
#endif
    }
    pfile.close();
    vfile.close();
    hfile.close();
#ifdef CATALOG
    hidfile.close();
#endif
  }

  //if (np > 0){
    cout<<" Ra min/max = "<<ramin<<"/"<<ramax<<endl;
    cout<<" Dec min/max = "<<decmin<<"/"<<decmax<<endl;
    cout<<" z min/max = "<<z_min<<"/"<<z_max<<endl;
    //}
  cout<<" file had "<<np<<" particles in volume."<<endl;

  return np;

}


//
// reads particles from file and returns vector of particle pointers
//
vector <Particle *> ReadMGSParticles(int &nread){
  nread = 0;

  cout<<"zmax: "<<zmax->GetVal()<<endl;

  vector <Particle*> particles;
  // Reserve max number of particles to avoid memory difficulties
  particles.reserve(MAXSIZE);
  cout<<"Reserved space for "<<MAXSIZE<<" particles."<<endl;

  cout<<"Reading MGS files."<<endl;

  for(int p=60;p<65;p++) {
    for(int t=60;t<65;t++) {
      if(particles.size()>MAXSIZE){
	cerr<<"[read_gadget] Exceeded maximum number of particles"
	    <<particles.size()<<" "<<MAXSIZE<<endl;	  exit(1);
      }
      cout<<"[read_MGS] Reading File "<<p<<" "<<t<<"..."<<endl;
      nread += ReadMGS(p,t,particles);
      cout<<nread<<endl;
    }
  }
  cout<<"Done."<<endl;
  PRNT("[read_MGS] end:",nread);
  assert(particles.size()>0);

  return particles;

}
