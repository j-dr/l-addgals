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
