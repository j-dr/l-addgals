//#include "cosmology.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "global_vars.h"

using namespace std;
const float amax = 0.1;

void Cosmology::GetZofR(float OmegaM, float OmegaL){
  int entries = 30000;
  //float amax = 0.1;
  double This_a, da, This_H;

  redshift.resize(entries);
  losdist.resize(entries);

  redshift[0] = 0;
  losdist[0] = 0;
  da = (1.0 - amax)/entries;

  for(int i=0;i<entries;i++)
    {
      This_a = 1.0 - da*i;
      This_H = 100.*sqrt(OmegaM/(This_a*This_a*This_a) + OmegaL);
      redshift[i] = 1.0/This_a - 1;
      losdist[i] = losdist[i-1] + 1.0/(This_H*This_a*This_a)*da*cspeed;
    }
  
}

double Cosmology::ZofR(double R){
  std::vector<double>::iterator pos;
  pos=upper_bound(losdist.begin(), losdist.end(), R);
  int i = distance(losdist.begin(),pos);
  double delta_x = losdist[i] - losdist[i-1];
  double dx      = R         - losdist[i-1];
  double delta_y = redshift[i] - redshift[i-1];
  double y = dx/delta_x*delta_y+redshift[i-1];

  return y;
}

double Cosmology::RofZ(double z){
  //if((z>1.4585)||(z<0)){
  if(z<0) z = 0; //changed by mbusha to work with Bolshoi snapshots
  //if((z>=4.0)||(z<0)){
  if((z>=1.0/amax - 1)||(z<0)){
    cerr<<"Redshift out of bounds in RofZ function: "<<z<<" "<<1.4585<<endl;
    exit(1);
  }
  std::vector<double>::iterator pos;
  pos=upper_bound(redshift.begin(), redshift.end(), z);
  int i = distance(redshift.begin(),pos);
  double delta_x = redshift[i] - redshift[i-1];
  double dx      = z        - redshift[i-1];
  double delta_y = losdist[i] - losdist[i-1];
  double y = dx/delta_x*delta_y+losdist[i-1];

  return y;
}



