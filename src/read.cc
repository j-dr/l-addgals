#include "particle.h"
#include "halo.h"
#include <vector>

vector <Halo*> ReadWHalos(void);
vector <Halo*> ReadHVHalos(void);
vector <Halo*> ReadGadgetHalos(void);
vector <Halo*> ReadRockstarHalos(void);
vector <Halo*> ReadRockstarHalosHlist(void);
vector <Halo*> ReadMGSHalos(void);
vector <Particle *> ReadWarren(int &nread,vector <Halo *> halos );
vector <Particle *> ReadHVParticles(int &nread);
vector <Particle *> ReadGadgetParticles(int &nread);
vector <Particle *> ReadMGSParticles(int &nread);


vector <Particle *> ReadParticles(int &nread, vector <Halo *> halos ){
  vector <Particle *> particles;

  //if(sim.Type() == "HV") particles = ReadHVParticles(nread);
  //else if(sim.Type() == "WAR") particles = ReadWarren(nread, halos);
  //else if(sim.Type() == "GADGET2") particles = ReadGadgetParticles(nread);
  if(sim.Type() == "GADGET2") particles = ReadGadgetParticles(nread);
  else if(sim.Type() == "MGS") particles = ReadMGSParticles(nread);
  else {cout<<"[ReadParticles arg arg] Wrong simulation type:"<<sim.Type()<<endl;exit(1);}
  return particles;
}



vector <Particle *> ReadParticles(int &nread){
  vector <Particle *> particles;
  nread = 0;
  /*
  if(sim.Type() == "HV") {
    particles = ReadHVParticles(nread);
  }
  */
  //else if(sim.Type() == "GADGET2"){
  if(sim.Type() == "GADGET2"){
    particles = ReadGadgetParticles(nread);
    //vector <Halo*> halos = ReadGadgetHalos();
  }
  /*
  else if(sim.Type() == "WAR"){
    vector <Halo*> halos = ReadWHalos();
    particles = ReadWarren(nread, halos);
  }
  */
  else if(sim.Type() == "MGS"){
    particles = ReadMGSParticles(nread);
    cout<<"Finished reading particles."<<endl;
  }
  else {
    cout<<"[ReadParticles] Wrong simulation type:"<<sim.Type()<<endl;
    exit(1);
  }
  return particles;
}


vector <Halo*> ReadHalos(void){
  vector <Halo *> halos;
  //if(sim.Type() == "HV") halos = ReadHVHalos();
  //else if(sim.Type() == "WAR") halos = ReadWHalos();
  //else if(sim.Type() == "GADGET2") halos = ReadGadgetHalos();
#ifdef SHAM_TEST
  if(sim.Type() == "GADGET2") halos = ReadRockstarHalosHlist();
#else
  if(sim.Type() == "GADGET2") halos = ReadRockstarHalos();
#endif
  else if(sim.Type() == "MGS") halos = ReadMGSHalos();
  else {cout<<"[ReadHalos] Wrong simulation type:"<<sim.Type()<<endl;exit(1);}
  return halos;
}



