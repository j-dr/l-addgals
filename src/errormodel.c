#include <time>
#include <math>
#include <random>
#include "errormodel.h"

using namespace std;

//Given simulated observed des magnitudes, apply an error model
//given exposure time, limiting mags, lnscatter constant over the sky 
//Returns:
//flux observed flux with errors included
//fluxerr observed flux errors
//omag observed magnitudes with errors included
//omagerr observed magnitude errors
void apply_uniform_errormodel(float exptime[], float limmags[], float lnscat, int nband, 
			      float zeropoint[], float nsigma, vector<float> &mag, 
			      vector<float> &flux, vector<float> &fluxerr, 
			      vector<float> &omag, vector<float> &omagerr)
{
  int i;
  float f1lim[nband];
  float fsky1[nband];
  float iexptime[nband];
  long int seed;
  seed = time(NULL);
  std::default_random_engine generator(seed);
  std::normal_distribution<float> rng();
  
  for (i=0; i<nband; i++)
    {
      f1lim[i] = pow(10, (limmags[i]-zeropoint[i])/(-2.5));
      fsky[i] = pow(f1lim[i],2)*exptimes[i]/pow(nsig,2) - f1lim[i];
      iexptime[i] = 1/exptime[i]
      if fsky[i] < 0.001 {
	  fsky[i] = 0.001;
	}
    }

  for (vector<float>::iterator itr=mag.begin(); itr!=mag.end; itr++)
    {
      i = itr-mag.begin();
      b = i%nband;
      flux[i] = exptimes[b] * pow( 10, 0.4 * ( *itr-zeropoint[b] ) );
      fluxerr[i] = sqrt( fsky1[b] * exptimes[b] + flux[i]);
      flux[i] += fluxerr[i] * rng(generator);
      fluxerr[i] += exp( log( fluxerr[i] + lnscat[b] * rng(generator) ) );
      flux[i] = flux[i] * iexptime[b];
      fluxerr[i] = fluxerr[i] * iexptime[b];
      omag[i] = zeropoint[b] - 2.5*log10(flux[i]/1.0e-9);
      omagerr[i] = 1.086*fluxerr[i]/flux[i];
    }
}

void add_des_photometric_errors(float maglim[], float sig, int nband,
				vector<float> &mag, vector<float> &flux,
				vector<float> &fluxerr, vector<float> &omag,
				vector<float> &omagerr)
{
  const float nk[] = {4.6, 12.9, 17.7, 45.1, 14.9}; //flux for sky at 24th magnitude
  float es[] = {9.95, 10.01, 8.02, 6.18, 0.812}; //flux for a 24th magnitude galaxy
  float exposure_time[nband];
  float zps[nband];
  seed = time(NULL);
  std::default_random_engine generator(seed);
  std::normal_distribution<float> rng();
  const float apperture = 1.5;   //angular size of a typical galaxy
  const float pixel = 0.27;      //angular size of a detector pixel
  const float apperture_area = M_PI * pow( apperture / 2.0, 2 );
  const float pixels_area = pow( pixel, 2 );
  const int npixels = apperture_area / pixels_area;


  //maglim = [24.0, 24.0, 24.0, 24.0, 24.0, 24.0]
  
  for (i=0; i<nband; i++)
    {
      zps[i] = 24. + 2.5 * log10( es[i] );  //zero-point for flux calculation of a galaxy
      es[i] = pow(10, 0.4 * ( zps[i] - maglim[i] ) );  //recalibrate our 10-sigma fluxes to detection limits
      exposure_time[i] = pow( 10.0 / es[i], 2 ) * ( nk[i] * npixels + es[i] );
    }

  for (vector<float>::iterator itr=mag.begin(); itr!=mag.end; itr++)
    {
      i = itr-mag.begin();
      b = i%nband;
      flux[i] = pow( 10.0, -0.4 * (*mag - zps[b]) ) * exposure_time[b];
      fluxerr[i] = sqrt( nk[b] * exposure_time[b] * npixels + flux[i] );
      flux[i] += rng(generator) * fluxerr[i];
      if (flux[i] < 0){
	omag[i] = 99.0;
      }	else {
	omag[i] = zps[b] - 2.5 * alog10( flux[i] / exposure_time[i] );
      }
      omagerr[i] = 2.5 * log10( 1.0 + fluxerr[i] / flux[i] );
    }
}
