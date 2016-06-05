#ifndef kcorrect_utils_h
#define kcorrect_utils_h
#include <cmath>
#include <istream>
#include <utility>
#include <iterator>
#include <vector>
#include <functional>
#include "kcorrect.h"

struct magtuple
{
  float bands[5];
};

struct ginfo
{
  int id, central;
  float refmag, redshift;
};

#ifdef UNITTESTS
std::istream & operator>>(std::istream & is, ginfo & in)
{
  float pass;
  is >> in.id >> in.refmag >> pass >> pass >>
    in.redshift >> pass;

  return is;
}
#endif

std::istream & operator>>(std::istream & is, magtuple & in);

template <class T> struct magop : std::binary_function <T,T,T> {
  T operator() (const T& lhs, const T& rhs) const {return 2.5*log10(lhs/rhs);}
};

template <class T> struct appmagnitude : std::binary_function <T,T,T> {
  T operator() (const T& maggie, const T& zeropoint) const {return zeropoint-2.5*log10(maggie/1.0e-9);}
};

template <class T> struct absmagnitude : std::binary_function <T,T,T> {
  T operator() (const T& appmag, const T& redshift) const {return appmag-z2dm(redshift,0.3,0.7);}
};

void reconstruct_maggies(float *coeff, float *redshift, int ngal, float zmin, 
			 float zmax, float band_shift, char filterfile[], float *maggies);

void match_coeff(std::vector<int> &sed_ids, float *coeffs);

void k_calculate_magnitudes(std::vector<float> &coeff, std::vector<float> &redshift,
			    float zmin, float zmax, float band_shift, int nband,
			    char filterfile[], std::vector<float> &omag,
			    std::vector<float> &amag);

void assign_colors(std::vector<float> &reference_mag, std::vector<float> &coeff, 
		   std::vector<float> &redshift, float zmin, float zmax,
		   float band_shift, int nbands, char filterfile[], 
		   std::vector<float> &omag, std::vector<float> &amag, 
		   std::vector<float> abcorr, bool refflag=true);


#endif
