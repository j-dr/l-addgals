#ifndef kcorrect_utils_h
#define kcorrect_utils_h
#include <cmath>
#include <istream>
#include <utility>
#include <iterator>
#include <vector>
#include <functional>

typedef std::pair<int,int> indpair;
bool paircomp (indpair l, indpair r)
{
return l.first<r.first; 
}

struct magtuple
{
  double bands[5];
};

struct ginfo
{
  int id;
  double refmag, redshift;
};

std::istream & operator>>(std::istream & is, magtuple & in)
{
  is >> in.bands[0] >> in.bands[1] >> in.bands[2] >>
    in.bands[3] >> in.bands[4];

  return is;
}

std::istream & operator>>(std::istream & is, ginfo & in)
{
  float pass;
  is >> in.id >> in.refmag >> pass >> pass >>
    in.redshift >> pass;

  return is;
}

template <class T> struct magnitude : std::binary_function <T,T,T> {
  T operator() (const T& lhs, const T& rhs) const {return 2.5*log10(lhs/rhs);}
};

template <class T> struct appmagnitude : std::binary_function <T,T,T> {
  T operator() (const T& maggie, const T& zeropoint) const {return zeropoint-2.5*log10(maggie/1.0e-9);}
};

void reconstruct_maggies(*float coeff, *float redshift, int ngal, float zmin, 
			 float zmax, float band_shift, char[] filterfile, *float maggies);

void match_coeff(vector<int> &sed_ids, double* coeffs);

void k_calculate_magnitudes(vector<double> &coeff, vector<double> &redshift,
			    float zmin, float zmax, float band_shift,
			    char[] filterfile, vector<double> &omag,
			    vector<double> &amag);

void assign_colors(vector<double> &reference_mag, vector<int> &sed_ids, 
		   vector<double> &redshift, float zmin, float zmax,
		   float band_shift, int nbands, char[] filterfile, 
		   vector<double> &omag, vector<double> &amag);

struct addgals_ginfo {
  
  int id, central;
  float
  
};

#endif
