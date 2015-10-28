#include "kcorrect_utils.h"
#include "kcorrect.h"
#include "cosmo.h"
#include <math.h>
#include <string>
#define FILESIZE 2000
#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}

using namespace std;

static int nz=1000;
static int nk,nv,nl,maxn;
static float *lambda=NULL;
static float *vmatrix=NULL;
static float *rmatrix=NULL;
static float *zvals=NULL;
static float *filter_lambda=NULL;
static float *filter_pass=NULL;
static int *filter_n=NULL;

void reconstruct_maggies(*float coeff, *float redshift, int ngal, float zmin, 
			 float zmax, float band_shift, char[] filterfile, *float maggies)
{

  int i,j,k,c,ndim,nchunk,ncurrchunk,*sizes=NULL;
  char vfile[FILESIZE],lfile[FILESIZE],ffile[FILESIZE],path[FILESIZE];
  char vmatrixfile[FILESIZE],lambdafile[FILESIZE],filterfile[FILESIZE];

  //Read in templates
  strcpy(vfile,"vmatrix.default.dat");
  strcpy(lfile,"lambda.default.dat");
  sprintf(path,"%s/data/templates",getenv("KCORRECT_DIR"));
  
  sprintf(vmatrixfile, "%s/%s", path, vfile);
  sprintf(lambdafile, "%s/%s", path, lfile);

  k_read_ascii_table(&vmatrix,&ndim,&sizes,vmatrixfile);
  nl=sizes[1];
  nv=sizes[0];
  FREEVEC(sizes);

  k_read_ascii_table(&lambda,&ndim,&sizes,lambdafile);
  if(sizes[0]!=nl+1) {
    fprintf(stderr,"vmatrix and lambda files incompatible (nl==%d vs sizes[0]=%d).\n",nl,sizes[0]);
    exit(1);
  }
  FREEVEC(sizes);

  //Load filters
  k_load_filters(&filter_n,&filter_lambda,&filter_pass,&maxn,&nk,filterfile);

  //Create matrix of projections of templates onto filters (rmatrix)
  rmatrix=(float *) malloc(nz*nv*nk*sizeof(float));
  zvals=(float *) malloc(nz*sizeof(float));
  for(i=0;i<nz;i++)
    zvals[i]=zmin+(zmax-zmin)*((float)i+0.5)/(float)nz;
  k_projection_table(rmatrix,nk,nv,vmatrix,lambda,nl,zvals,nz,filter_n,
		     filter_lambda,filter_pass,band_shift,maxn);

  k_reconstruct_maggies(zvals,nz,rmatrix,nk,nv,*coeff,*redshift,*maggies,ngal);

}

void match_coeff(vector<int> &sed_ids, double* coeffs)
{

  static int ntemps=5;
  int n;
  string filename = "../training/cooper/combined_dr6_cooperdat";
  istream coeff_file(filename.c_str());
  if (coeff_file.fail()) {
    cerr<<"error: cannot open "<<filename<<endl;
    exit(1);
  }

  vector<magtuple> sdsscoeffs(sed_ids.size());
  
  copy(istream_iterator<magtuple>(coeff_file),
	    istream_iterator<magtuple>(),
	    back_inserter(sdsscoeffs));
    
  for (vector<int>::iterator it = sed_ids.begin;
       it!=sed_ids.end(), it++)
    {
      for (n=0; n<ntemps, ++n) {
	coeffs[it*ntemps+n] = sdsscoeffs[*it].coeffs[n];
      }
    }
}

void k_calculate_magnitudes(vector<double> &coeff, vector<double> &redshift,
			    float zmin, float zmax, float band_shift,
			    char[] filterfile, vector<double> &omag,
			    vector<double> &amag)
{
  int ngal = redshift.size();
  vector<double> kcorrection(ngal);

  //observer frame maggies
  reconstruct_maggies(&coeff[0], &redshift[0], ngal, zmin,
                      zmax, 0.0, filterfile, &omag[0]);

  //rest frame maggies
  reconstruct_maggies(&coeff[0], &redshift[0], ngal, zmin,
                      zmax, band_shift, filterfile, &omag[0]);

  //calculate k-correction in mags
  transform(omag.begin(), omag.end(), amag.begin(), 
		 kcorrection.begin(), magnitude());

  //calculate apparent magnitudes
  transform(omag.begin(), omag.end(), omag.begin(),
		 appmagnitude());

  //calculate absolute magnitudes without kcorrection
  transform(omag.begin(), omag,end(), redshift.begin(),
		 amag.begin(), absmag());

  //add k-correction to absolute magnitudes
  transform(amag.begin(), amag.end(), kcorrection.begin(),
		 amag.begin(), plus());

}

void assign_colors(vector<double> &reference_mag, vector<int> &sed_ids, 
		   vector<double> &redshift, float zmin, float zmax,
		   float band_shift, int nbands, char[] filterfile, 
		   vector<double> &omag, vector<double> &amag){
  int ngal = reference_mag.size();
  int ntemp = 5;
  double sdss_bandshift = 0.1;
  char sdss_filterfile[FILESIZE];
  vector<double> coeff(ngal*ntemp);
  vector<double> deltam(ngal);

  strcopy(sdss_filterfile, "./sdss_filter.txt");

  //Match the correct coefficients to each galaxy
  match_coeff(sed_ids, &coeff[0]);

  //calculate SDSS r-band abs mags to determine magnitude shifts to apply
  //to kcorrect outputs for other bands
  k_calculate_magnitudes(&coeff[0], &redshift[0], zmin, zmax, sdss_bandshift,
			 sdss_filterfile, &omag[0], &amag[0]);

  //determine shift needed to get reference_mag from sdss_amag
  transform(reference_mag.begin(), reference_mag.end(), amag.begin(),
	    deltam.begin(), minus());
    
  //calculate observed and abs mags in desired output bands without shift
  //calculated from SDSS bands
  k_calculate_magnitudes(&coeff[0], &redshift[0], zmin, zmax, band_shift,
			 filterfile, &omag[0], &amag[0]);
  
  //apply SDSS shift to app and abs mags (same for all bands)
  transform(omag.begin(), omag.end(), deltam.begin(), omag.begin(), 
	    plus());
  transform(amag.begin(), amag.end(), deltam.begin(), amag.begin(), 
	    plus());

}

int main(int argc, char **argv)
{
  if (argc < 2){
    cerr << "Usage: " << argv[0] << " ginfo1.dat" << endl;
    return 1;
  }

  istream ginfo_file(argv[1]);
  float zmin, zmax, band_shift;
  int nbands;
  char filter_file[FILESIZE];
  vector<ginfo> sdss_ginfo;

  if (ginfo_file.fail()) {
    cerr<<"error: cannot open "<<argv[1]<<endl;
    exit(1);
  }
  
  zmin = 0.0;
  zmax = 0.5;
  band_shift = 0.1;
  nbands = 5;
  strcpy(filter_file, "./des_filters.txt");

  copy(istream_iterator<ginfo>(ginfo_file),
	    istream_iterator<ginfo>(),
	    back_inserter(sdss_ginfo));

  int ngal = sdss_ginfo.size();
  vector<double> sed_ids(ngal);
  vector<double> redshift(ngal);
  vector<double> refmag(ngal);
  vector<double> omag(ngal*nbands);
  vector<double> amag(ngal*nbands);

  typedef vector<ginfo>::iterator gitr;
  for (gitr it sdss_ginfo.begin(), it!=sdss_ginfo.end();++it)
    {
      sed_ids[it] = it->id;
      redshift[it] = it->redshift;
      refmag[it] = it->refmag;
    }

  assign_colors(&refmag[0], &sed_ids[0], &redshift[0], zmin, zmax,
		band_shift, nbands, filterfile, &omag, &amag);

  //write out magnitudes
  ofstream omag_file("./des_omagtest.txt");
  ofstream amag_file("./des_amagtest.txt");
  typedef vector<double>::iterator ditr;

  for (ditr it=omag.begin(); it!=omag.end(), ++it)
    {
      omag_file << *it;
      if ((it!=omag.begin()) && (it%nbands==0)) omag_file << "\n";
    }

    for (it=amag.begin(); it!=amag.end(), ++it)
    {
      amag_file << *it;
      if ((it!=amag.begin()) && (it%nbands==0)) amag_file << "\n";
    }

}
