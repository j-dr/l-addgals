#include "hv.h"
#include "hv2.h"
#include "functions.h"
#include <fstream>
#define COLORS
#define PRINTHALOS
#include "chainel.h"
//#define NOCFS


#ifdef HVL
Simulation sim(box, use_cells);
#endif
#ifdef LANL
Simulation sim(box);
#endif

vector <Particle*> ReadParticles(int &);
Cosmology cosmo = sim.SimCosmology();
using namespace std;

int main(void){
  cosmo.Print();
  cout<<sim.Boxsize()<<endl;
  cout<<sim.Boxsize()*BOXFR<<endl;
  cosmo.ReadZFile();
  int nread;
  vector <Particle*> particles = ReadParticles(nread);
  //  float volume_fraction = particles.size()*1.0/nread;
  cout<<" Read "<<particles.size()<<" particles"<<endl;
  assert(particles.size()>0);
  ofstream outfile("hv_dm.dat");
  for(int i=0;i<particles.size();i++){
    if(randbool(0.5))
      particles[i]->Pos2Write(outfile);
  }
}
