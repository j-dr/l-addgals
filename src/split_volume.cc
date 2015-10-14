#include "box.h"
#include "constants.h"
#include "cosmo.h"
#include <sstream>
#include <string>

#ifdef PARALLEL
//void split_volume(int NTasks, int ThisTask)
void split_volume(void)
{
#define nzbins 2000

  float r_arr[nzbins];
  float z_arr[nzbins];
  float volume_arr[nzbins];
  float del_z;
  float PI = 3.1415926535897932384626433832795;
  float decmin_rad = PI*(90.0 - DECMAX)/180.0;
  float decmax_rad = PI*(90.0 - DECMIN)/180.0;
  float ramin_rad = PI*RAMIN/180.0;
  float ramax_rad = PI*RAMAX/180.0;
  float dra = ramax_rad - ramin_rad;
  float r0 = cosmo.RofZ(SIM_ZREDMIN);

  del_z = (SIM_ZREDMAX - SIM_ZREDMIN) / float(nzbins);

  for(int i=0;i<nzbins;i++){
    z_arr[i] = SIM_ZREDMIN + i*del_z;
    r_arr[i] = cosmo.RofZ(z_arr[i]);
    volume_arr[i] = (ramax_rad-ramin_rad)*(cos(decmin_rad)-cos(decmax_rad))*(pow(r_arr[i],3)-pow(r0,3))/3.;
  }
  float del_v = volume_arr[nzbins-1]/float(NTasks);
  float This_vmin = del_v*ThisTask;
  float This_vmax = This_vmin + del_v;

  ZREDMIN = 0.0;
  ZREDMAX = 0.0;
  for (int i=0;i<nzbins;i++)
    if (volume_arr[i] > This_vmin){
      ZREDMIN = z_arr[i];
      if (i > 0) 
	ZREDMIN = z_arr[i-1];
      break;
    }
  for (int i=0;i<nzbins;i++){
    if (volume_arr[i] > This_vmax){
      ZREDMAX = z_arr[i];
      break;
    }
    if (i == nzbins-1)
      ZREDMAX = z_arr[i];
  }
  //Magmin = -13.1 - 8.66*ZREDMIN + 2.9*ZREDMIN*ZREDMIN; //calculate our interpolating function
  Magmin = -14.2 - 22.0*ZREDMIN + 28.*ZREDMIN*ZREDMIN;
  cout<<"Task "<<ThisTask<<" is looking at the region z = "<<ZREDMIN<<"-"<<ZREDMAX<<endl;
}
#endif
