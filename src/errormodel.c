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
//flux    -- observed flux with errors included
//fluxerr -- observed flux errors
//omag    -- observed magnitudes with errors included
//omagerr -- observed magnitude errors
void apply_uniform_errormodel(float exptime[], float limmags[], float lnscat[], 
			      int nband, float zeropoint[], float nsigma, 
			      vector<float> &mag, vector<float> &flux, 
			      vector<float> &ivar, vector<float> &omag, 
			      vector<float> &omagerr)
{
  int i,b;
  float f1lim[nband];
  float fsky1[nband];
  float iexptime[nband];
  gsl_rng *rng;
  //long int seed;
  unsigned int seed;
  FILE *urandom;

  urandom = fopen ("/dev/urandom", "r");
  fread (&seed, sizeof (seed), 1, urandom);
  //seed = (long int)time(NULL);
  rng = gsl_rng_alloc(gsl_rng_ranlxd2);
  gsl_rng_set(rng, seed);
  
  for (i=0; i<nband; i++)
    {
      f1lim[i] = pow(10.0, (limmags[i]-zeropoint[i])/(-2.5));
      if (f1lim[i] < 120.0/exptime[i]) {
	  f1lim[i] = 120.0/exptime[i];
      }
      fsky1[i] = pow(f1lim[i],2)*exptime[i]/pow(nsigma,2) - f1lim[i];
      iexptime[i] = 1/exptime[i];
    }

  for (vector<float>::iterator itr=mag.begin(); itr!=mag.end(); itr++)
    {
      i = itr-mag.begin();
      b = i%nband;
      flux[i] = exptime[b] * pow( 10, -0.4 * ( *itr-zeropoint[b] ) );
      ivar[i] = exp( log( sqrt( fsky1[b] * exptime[b] + flux[i] ) )
			+ lnscat[b] * gsl_ran_gaussian(rng, 1.0) );
      flux[i] += ivar[i] * gsl_ran_gaussian(rng, 1.0);
      flux[i] = flux[i] * iexptime[b];
      ivar[i] = ivar[i] * iexptime[b];
      omag[i] = zeropoint[b] - 2.5 * log10( flux[i] );
      omagerr[i] = 1.086*ivar[i]/flux[i];
      if ((omag[i]>99.0) || (omag[i]!=omag[i])) {
	omag[i] = 99.0;
	omagerr[i] = 99.0;
      }
      ivar[i] = 1/(ivar[i]*ivar[i]);
    }
}

void observe_des_y5(vector<float> &mag, vector<float> &flux, 
		    vector<float> &ivar, vector<float> &omag,
		    vector<float> &omagerr, vector<bool> &idx)
{
  float maglim[] = {24.956,24.453,23.751,23.249,21.459};
  float maglim_cut[] = {25.5, 25.0, 24.4, 23.9, 22.0};
  float zeropoint[] = {22.5, 22.5, 22.5, 22.5, 22.5};
  float exptime[] = {14467.00,12471.00,6296.00,5362.00,728.00};
  float lnscat[] = {0.2,0.2,0.2,0.2,0.2};
  float delta_maglim = 2.0;
  int nband = 5;
  int ngal = mag.size()/nband;
  int count = 0;
  bool good=false;
  int i, b;

  for (vector<float>::iterator itr=mag.begin(); itr!=mag.end(); itr++)
    {
      i = (itr-mag.begin())/nband;
      b = (itr-mag.begin())%nband;
      good = good || (*itr <= (maglim_cut[b] + delta_maglim));
      if (b==(nband-1)) {
	if (good) {
	  copy(itr-nband+1, itr+1, mag.begin()+count);
	  count+=nband;
	}
	idx[i] = good;
	good = false;
      }
    }

  mag.resize(count);
  flux.resize(count);
  ivar.resize(count);
  omag.resize(count);
  omagerr.resize(count);

  apply_uniform_errormodel(exptime, maglim, lnscat, nband, zeropoint,
			   10.0, mag, flux, ivar, omag, omagerr);

}

#ifndef BCC
int main(int argc, char *argv[])
{
  int i, ngal;
  int nband = 5;
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
  vector<float> omag(ngal*nband,-99);
  vector<float> omagerr(ngal*nband);
  vector<bool> idx(ngal);
  
  observe_des_y5(mag, flux, fluxerr, omag, omagerr, idx);
  
  ofstream error_file("./des_errortest.txt");
  ofstream magcut_file("./des_magcut.txt");
  vector<float>::iterator ditr;

  for (ditr=omag.begin(); ditr!=omag.end(); ++ditr)
    {
      i = ditr-omag.begin();
      error_file << flux[i] << " " << fluxerr[i] << " " 
		 << *ditr << " " << omagerr[i] << endl;

      if ( (i+1) % nband == 0 )
	{
	  error_file << "\n";
	} else error_file << " ";
    }

  for (ditr=mag.begin(); ditr!=mag.end(); ++ditr)
    {
      i = ditr-mag.begin();
      magcut_file << *ditr;

      if ( (i+1) % nband == 0 )
	{
	  magcut_file << "\n";
	} else magcut_file << " ";
    }
}
#endif
