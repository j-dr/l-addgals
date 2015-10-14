#include "simulation.h"
#include "box.h"
#include <string>
#include <vector>

extern Simulation sim;
using namespace std;

/*
double LengthUnit(){
  double lengthunit;
  if (sim.Type()==HV)
    lengthunit=3000.0;
  else
    lengthunit=sim.Boxsize();
  return lengthunit;
}
*/

string Dir(){
  string dir;
  if(sim.Type()==WAR){
    dir = datadir+"warren/";
  }
  else if(sim.Type()==HV){
#ifdef HVL
    if(box == HVLIGHT){
      dir = datadir+"hubble/";
    }
    else if(box == HVZ0){
     dir = datadir+"hubble/z0/";
    }
    else if(box == HVPOW){
      //      dir= "/data2/hubble/POW/rnn/";
      // dir= datadir+"hubble/POW/";
      // dir= datadir+"hubble/POW2/";
       dir= datadir+"hubble/PO/";
    }
    else{
      cerr<<"[error:] Inconsistent box/sim type"<<simulation<<"\t"<<box<<endl;
      exit(1);
    }
#endif
  }
  else if(sim.Type()==GADGET2){
    dir = datadir;
  }
  else if (sim.Type()==MGS){
    dir = datadir;
  }
  else{
    cerr<<"[error:] Inconsistent box/sim type"<<simulation<<"\t"<<box<<endl;
    exit(1);
  }
  return dir;
}

/*double area(double alpha, double* dummy){
  return cos(alpha*3.141592654/180.);
}

double fractional_area(){
  
  integrator area_integrator(&area, 0,1.0e-03,0.01);
  double sky_coverage = (DECMAX-DECMIN)*area_integrator.Integrate(RAMIN,RAMAX);
  return sky_coverage/5156.62;
}
*/

