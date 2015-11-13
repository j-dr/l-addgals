#include <iostream>
#include <fstream>
#include <iterator>
#include <ctime>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "errormodel.h"

using namespace std;

//Given simulated observed des magnitudes, apply an error model
//given exposure time, limiting mags, lnscatter constant over the sky 
//Returns:
//flux observed flux with errors included
//fluxerr observed flux errors
//omag observed magnitudes with errors included
//omagerr observed magnitude errors
void apply_uniform_errormodel(float exptime[], float limmags[], float lnscat[], int nband, 
			      float zeropoint[], float nsigma, vector<float> &mag, 
			      vector<float> &flux, vector<float> &fluxerr, 
			      vector<float> &omag, vector<float> &omagerr)
{
  int i,b;
  float f1lim[nband];
  float fsky1[nband];
  float iexptime[nband];
  gsl_rng *rng;
  long int seed;
  seed = (long int)time(NULL);
  rng = gsl_rng_alloc(gsl_rng_ranlxd2);
  gsl_rng_set(rng, seed);
  
  for (i=0; i<nband; i++)
    {
      f1lim[i] = pow(10, (limmags[i]-zeropoint[i])/(-2.5));
      fsky1[i] = pow(f1lim[i],2)*exptime[i]/pow(nsigma,2) - f1lim[i];
      iexptime[i] = 1/exptime[i];
      if (fsky1[i] < 0.001) {
	  fsky1[i] = 0.001;
	}
    }

  for (vector<float>::iterator itr=mag.begin(); itr!=mag.end(); itr++)
    {
      i = itr-mag.begin();
      b = i%nband;
      flux[i] = exptime[b] * pow( 10, 0.4 * ( *itr-zeropoint[b] ) );
      fluxerr[i] = sqrt( fsky1[b] * exptime[b] + flux[i]);
      flux[i] += fluxerr[i] * gsl_ran_gaussian(rng, 1.0);
      fluxerr[i] += exp( log( fluxerr[i] + lnscat[b] * gsl_ran_gaussian(rng, 1.0) ) );
      flux[i] = flux[i] * iexptime[b];
      fluxerr[i] = fluxerr[i] * iexptime[b];
      omag[i] = zeropoint[b] - 2.5*log10(flux[i]/1.0e-9);
      omagerr[i] = 1.086*fluxerr[i]/flux[i];
    }
}

void add_des_photometric_errors(float maglim[], int nband, vector<float> &mag,
				vector<float> &flux, vector<float> &fluxerr, 
				vector<float> &omag, vector<float> &omagerr,
				vector<bool> &good)
{
  const float nk[] = {4.6, 12.9, 17.7, 45.1, 14.9}; //flux for sky at 24th magnitude
  float es[] = {9.95, 10.01, 8.02, 6.18, 0.812}; //flux for a 24th magnitude galaxy
  float exposure_time[nband];
  float zps[nband];
  gsl_rng *rng;
  long int seed;
  seed = (long int)time(NULL);
  rng = gsl_rng_alloc(gsl_rng_ranlxd2);
  gsl_rng_set(rng, seed);
  const float apperture = 1.5;   //angular size of a typical galaxy
  const float pixel = 0.27;      //angular size of a detector pixel
  const float apperture_area = M_PI * pow( apperture / 2.0, 2 );
  const float pixels_area = pow( pixel, 2 );
  const int npixels = apperture_area / pixels_area;
  int i,b;

  //maglim = [24.0, 24.0, 24.0, 24.0, 24.0, 24.0]
  
  for (i=0; i<nband; i++)
    {
      zps[i] = 24. + 2.5 * log10( es[i] );  //zero-point for flux calculation of a galaxy
      es[i] = pow(10, 0.4 * ( zps[i] - maglim[i] ) );  //recalibrate our 10-sigma fluxes to detection limits
      exposure_time[i] = pow( 10.0 / es[i], 2 ) * ( nk[i] * npixels + es[i] );
    }

  for (vector<float>::iterator itr=mag.begin(); itr!=mag.end(); itr++)
    {
      i = itr-mag.begin();
      b = i%nband;
      flux[i] = pow( 10.0, -0.4 * (*itr - zps[b]) ) * exposure_time[b];
      fluxerr[i] = sqrt( nk[b] * exposure_time[b] * npixels + flux[i] );
      flux[i] += gsl_ran_gaussian(rng, 1.0) * fluxerr[i];
      if (flux[i] < 0){
	omag[i] = 99.0;
      }	else {
	omag[i] = zps[b] - 2.5 * log10( flux[i] / exposure_time[i] );
      }
      omagerr[i] = 2.5 * log10( 1.0 + fluxerr[i] / flux[i] );
      good[i] =  good[i] | (omag[i] > maglim[b]) ? true : false;
    }
}

int main(int argc, char *argv[])
{
  int i, ngal;
  int nband = 5;
  float maglim[] = {24.0, 24.0, 24.0, 24.0, 24.0};
  ifstream mag_file(argv[1]);
  vector<float> mag;
  
  if (mag_file.fail())
    {
      cout << "error: cannot open file" << endl;
    }

  copy(istream_iterator<float>(mag_file),
	    istream_iterator<float>(),
	    back_inserter(mag));  

  ngal = mag.size() / nband;
  vector<float> flux(ngal*nband);
  vector<float> fluxerr(ngal*nband);
  vector<float> omag(ngal*nband);
  vector<float> omagerr(ngal*nband);
  vector<bool> good(ngal, false);
  
  add_des_photometric_errors(maglim, nband, mag, flux, fluxerr, 
			     omag, omagerr, good);

  ofstream error_file("./des_errortest.txt");
  vector<float>::iterator ditr;

  for (ditr=omag.begin(); ditr!=omag.end(); ++ditr)
    {
      i = ditr - omag.begin();
      if ( !good[i/nband] ) continue;

      error_file << flux[i] << " " << fluxerr[i] << " " 
		 << *ditr << " " << omagerr[i] << endl;

      if ( (i+1) % nband == 0 )
	{
	  error_file << "\n";
	} else error_file << " ";
    }
}
