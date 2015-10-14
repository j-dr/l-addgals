#include <vector>
#include <cassert> 
#include <unistd.h>     //for sleep()
#include <fstream>
#include "hv2.h"
#include "hv.h"

#define COLORS
#define PRINTHALOS
//#define NOCFS


#ifdef HVL
Simulation sim(box, use_cells);
#endif
#ifdef LANL
Simulation sim(box);
#endif
#ifdef GADGET
Simulation sim(box);
#endif

Cosmology cosmo = sim.SimCosmology();
using namespace std;

extern "C" void covar(float *x1, float *y1, float *z1, 
		      float *vx1, float *vy1, float *vz1, int np1,
		      float rcube, float rmin, float rmax, int nbin);
extern double normal_random(float mean, float stddev);
extern double ranf(void);


//extern vector <Galaxy *> GetGalaxies(double vol, ChainEl chel);

int main(void){
  cosmo.Print();
  vector <Galaxy *> galaxies;
  ChainEl chainel(-0.475, -0.2517, 2.856, 0.0, 2.2, -0.0536, 0.8785); 

  float volume = sim.LengthUnit()*sim.LengthUnit()*sim.LengthUnit();
  cout<<"Volume in lightcone: "<<volume<<endl;
  galaxies=GetGalaxies(volume, chainel);
}
