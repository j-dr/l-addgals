#include <vector>
#include "outputs.h"
#include "kcorrect_utils.h"
#include "cosmo.h"
#include <math.h>
#include <cstring>
#include <limits>
#include <string>
#include <iomanip>
#include <algorithm>

Cosmology cosmo(0.3, 0.82, 0.7);
#define FILESIZE 2000
using namespace std;

bool ccfits_test()
 {
   
   int size = 1000;
   vector<Galaxy *> galaxies(size);
   vector<Particle *> particles(size);
   vector<float> amag(size*5);
   vector<float> tmag(size*5);
   vector<float> mr(size);
   vector<float> omag(size*5);
   vector<float> omagerr(size*5);
   vector<float> flux(size*5);
   vector<float> ivar(size*5);
   vector<double> e(size*2);
   vector<double> s(size);
   vector<bool> idx(size,true);
   vector<Halo *> halos(size);
   vector<int> sed_ids(size);
   vector<float> coeffs(size*5);
   string outgfn = "test.fits";
   string outghfn = "testh.fits";
   bool complete = false;
   int i;
   Point pt;
   vector<Particle *> p;
   Particle ptc(pt, pt, 1.0);
   p.push_back(&ptc);
   Galaxy g(1.0);
   g.P(p[0]);

   cout << "filling in gals" <<endl;
   for (i=0; i<size; i++) {
     galaxies[i] = &g;
   }

   cout << "Writing test catalog to fits" << endl;
   try {
     write_bcc_catalogs(galaxies, particles, amag, tmag, mr,
			omag, omagerr, flux, ivar, e, s,
			idx, halos, sed_ids, coeffs,
			outgfn, outghfn);
     complete = true;
   }
   catch (...) {
     complete = false;
   }

   return complete;
}

int kcorrect_test()
{
  string ginfo_name = "/nfs/slac/g/ki/ki21/cosmo/jderose/addgals/catalogs/Buzzard/Lb1050_v1.1/0/000/hv_output/gal_ginfo1.dat";
  cosmo.GetZofR(cosmo.OmegaM(),cosmo.OmegaL());
  ifstream ginfo_file(ginfo_name.c_str());
  float zmin, zmax, band_shift;
  int nbands, ntemp, i;
  char filterfile[FILESIZE];
  vector<ginfo> sdss_ginfo;
  bool complete = false;

  if (ginfo_file.fail()) {
    cerr<<"error: cannot open "<<ginfo_name<<endl;
    exit(1);
  }
  
  zmin = 0.0;
  //zmax = 0.5;
  zmax=2.0;
  band_shift = 0.1;
  nbands = 6;
  strcpy(filterfile, "./des_filters.txt");

  copy(istream_iterator<ginfo>(ginfo_file),
	    istream_iterator<ginfo>(),
	    back_inserter(sdss_ginfo));

  int ngal = sdss_ginfo.size();
  ntemp = 5;
  vector<int> sed_ids(ngal);
  vector<float> redshift(ngal);
  vector<float> refmag(ngal);
  vector<float> omag(ngal*nbands);
  vector<float> amag(ngal*nbands);
  vector<float> coeff(ngal*ntemp);
  vector<ginfo>::iterator itr;

  for (itr=sdss_ginfo.begin();itr<sdss_ginfo.end();++itr)
    {
      sed_ids[itr-sdss_ginfo.begin()] = itr->id;
      redshift[itr-sdss_ginfo.begin()] = itr->redshift;
      refmag[itr-sdss_ginfo.begin()] = itr->refmag;
    }
  cout<<"First few sdss sed ids: "<<endl;
  for (i=0; i<5; i++){
    cout<<sed_ids[i]<<", ";
  }
  cout<<endl;


  match_coeff(sed_ids, &coeff[0]);
  assign_colors(refmag, coeff, redshift, zmin, zmax,
		band_shift, nbands, filterfile, omag, amag);

  //write out magnitudes
  ofstream omag_file("./des_omagtest.txt");
  ofstream amag_file("./des_amagtest.txt");
  vector<float>::iterator ditr;

  for (ditr=omag.begin(); ditr!=omag.end(); ++ditr)
    {
      omag_file << *ditr;
      if ((ditr-omag.begin()+1)%nbands==0)
	{
	  omag_file << "\n";
	} else omag_file << " ";
    }

  for (ditr=amag.begin(); ditr!=amag.end(); ++ditr)
    {
      amag_file << *ditr;
      if ((ditr-amag.begin()+1)%nbands==0)
	{
	  amag_file << "\n";
	} else amag_file << " ";
    }
  
  complete = true;
  return complete;

}

int main(int argc, char **argv)
{

  if (!ccfits_test()) {
    cout << "ccfits test: FAILED!" << endl;
  }
  else {
    cout << "ccfits test: PASSED!" <<endl;
  }

  if (!kcorrect_test()) {
    cout << "kcorrect test: FAILED!" << endl;
  }
  else {
    cout << "kcorrect test: PASSED!" <<endl;
  }

}




