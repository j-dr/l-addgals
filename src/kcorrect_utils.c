#include "kcorrect_utils.h"
#include "kcorrect.h"
#include "cosmo.h"
#include <math.h>
#include <cstring>
#include <limits>
#include <string>
#include <iomanip>
#include <algorithm>
#define FILESIZE 2000
#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}

Cosmology cosmo(0.286, 0.82, 0.7);
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

void reconstruct_maggies(float *coeff, float *redshift, int ngal, float zmin, 
			 float zmax, float band_shift, char filterfile[], float *maggies)
{
  cout<<"Reconstructing maggies"<<endl;

  int i,ndim,*sizes=NULL;
  char vfile[FILESIZE],lfile[FILESIZE],path[FILESIZE];
  char vmatrixfile[FILESIZE],lambdafile[FILESIZE];

  cout<<"First few coefficients: "<<endl;
  for (i=0;i<5;i++){
    cout<<coeff[i]<<", ";
  }
  cout<<endl;

  //Read in templates
  strcpy(vfile,"vmatrix.default.dat");
  strcpy(lfile,"lambda.default.dat");
  sprintf(path,"%s/data/templates",getenv("KCORRECT_DIR"));
  
  sprintf(vmatrixfile, "%s/%s", path, vfile);
  sprintf(lambdafile, "%s/%s", path, lfile);

  cout<<"Reading in vmatrix"<<endl;
  k_read_ascii_table(&vmatrix,&ndim,&sizes,vmatrixfile);
  nl=sizes[1];
  nv=sizes[0];
  FREEVEC(sizes);

  cout<<"Reading in lambda file"<<endl;
  k_read_ascii_table(&lambda,&ndim,&sizes,lambdafile);
  if(sizes[0]!=nl+1) {
    fprintf(stderr,"vmatrix and lambda files incompatible (nl==%d vs sizes[0]=%d).\n",nl,sizes[0]);
    exit(1);
  }
  FREEVEC(sizes);

  //Load filters
  cout<<"Loading filters from "<<filterfile<<endl;
  k_load_filters(&filter_n,&filter_lambda,&filter_pass,&maxn,&nk,filterfile);
  cout<<"Allocating memory for projection table"<<endl;

  //Create matrix of projections of templates onto filters (rmatrix)
  rmatrix=(float *) malloc(nz*nv*nk*sizeof(float));
  zvals=(float *) malloc(nz*sizeof(float));
  cout<<"Creating projection table"<<endl;
  for(i=0;i<nz;i++)
    zvals[i]=zmin+(zmax-zmin)*((float)i+0.5)/(float)nz;
  k_projection_table(rmatrix,nk,nv,vmatrix,lambda,nl,zvals,nz,filter_n,
		     filter_lambda,filter_pass,band_shift,maxn);

  cout<<"Calling k_reconstruct_maggies"<<endl;
  k_reconstruct_maggies(zvals,nz,rmatrix,nk,nv,coeff,redshift,maggies,ngal);
  cout<<"Done reconstructing maggies"<<endl;
  cout<<"First few maggies: "<<endl;
  for (i=0;i<5;i++){
    cout<<maggies[i]<<", ";
  }
  cout<<endl;

}

void match_coeff(vector<int> &sed_ids, double* coeffs)
{
  cout<<"Matching coeffs"<<endl;
  static int ntemps=5;
  int n;
  string filename = "../training/cooper/combined_dr6_cooper.dat";
  ifstream coeff_file(filename.c_str());
  if (coeff_file.fail()) {
    cerr<<"error: cannot open "<<filename<<endl;
    exit(1);
  }

  vector<magtuple> sdsscoeffs;
  vector<int>::iterator it;

  copy(istream_iterator<magtuple>(coeff_file),
	    istream_iterator<magtuple>(),
	    back_inserter(sdsscoeffs));
  cout.precision(15);
  cout<<"First few sdss coefficients: "<<endl;
  for (n=0;n<ntemps;n++){
    cout<<setprecision(numeric_limits<double>::digits10+1)
	<<sdsscoeffs[0].bands[n]<<", ";
  }
  cout<<endl;

  
  for (it = sed_ids.begin(); it!=sed_ids.end(); it++)
    {
      for (n=0; n<ntemps; ++n) {
	coeffs[(it-sed_ids.begin())*ntemps+n] = sdsscoeffs[*it].bands[n];
      }
    }
  cout<<"Done matching coeffs"<<endl;
  cout<<"First few matches: "<<endl;
  for (n=0;n<5;n++){
    cout<<coeffs[n]<<", ";
  }
  cout<<endl;

}

void k_calculate_magnitudes(vector<double> &coeff, vector<double> &redshift,
			    float zmin, float zmax, float band_shift,
			    char filterfile[], vector<double> &omag,
			    vector<double> &amag)
{
  cout<<"Calculating magnitudes"<<endl;
  int i;
  int ngal = redshift.size();
  int nband = omag.size()/ngal;
  cout<<"Number of bands: "<<nband<<endl;
  vector<double> kcorrection(omag.size());

  //observer frame maggies
  reconstruct_maggies((float*)&coeff[0], (float*)&redshift[0], ngal, zmin,
                      zmax, 0.0, filterfile, (float*)&omag[0]);

  //rest frame maggies
  reconstruct_maggies((float*)&coeff[0], (float*)&redshift[0], ngal, zmin,
                      zmax, band_shift, filterfile, (float*)&amag[0]);

  cout<<"Calculating kcorrections"<<endl;
  //calculate k-correction in mags
  transform(omag.begin(), omag.end(), amag.begin(), 
	    kcorrection.begin(), magop<double>());

  cout<<"Calculating apparent magnitudes"<<endl;
  //calculate apparent magnitudes
  transform(omag.begin(), omag.end(), omag.begin(),
	    bind1st(appmagnitude<double>(),22.5));

  cout<<"Calculating absolute magnitudes"<<endl;
  //calculate absolute magnitudes without kcorrection

  vector<double>::iterator it;
  for (it=redshift.begin(); it!=redshift.end(); it++){
    for (i=0; i<nband; i++){
      amag[nband*(it-redshift.begin())+i] = AbsMag(omag[nband*(it-redshift.begin())+i],*it);
    }
  }

	//transform(omag.begin(), omag.end(), redshift.begin(),
	//	    amag.begin(), AbsMag);

  cout<<"Adding kcorrection"<<endl;
  //add k-correction to absolute magnitudes
  transform(amag.begin(), amag.end(), kcorrection.begin(),
	    amag.begin(), plus<double>());

  cout<<"Done calculating magnitudes"<<endl;

}

void assign_colors(vector<double> &reference_mag, vector<int> &sed_ids, 
		   vector<double> &redshift, float zmin, float zmax,
		   float band_shift, int nbands, char filterfile[], 
		   vector<double> &omag, vector<double> &amag){
  int i;
  int ngal = reference_mag.size();
  int ntemp = 5;
  double sdss_bandshift = 0.1;
  char sdss_filterfile[FILESIZE];
  vector<double> coeff(ngal*ntemp);
  vector<double> deltam(ngal);

  //sprintf(sdss_filterfile,"%s/data/templates/sdss_filters.dat",getenv("KCORRECT_DIR"));
  strcpy(sdss_filterfile, "./sdss_filter.txt");

  //Match the correct coefficients to each galaxy
  match_coeff(sed_ids, &coeff[0]);

  //calculate SDSS r-band abs mags to determine magnitude shifts to apply
  //to kcorrect outputs for other bands
  k_calculate_magnitudes(coeff, redshift, zmin, zmax, sdss_bandshift,
			 sdss_filterfile, omag, amag);

  //determine shift needed to get reference_mag from sdss_amag
  transform(reference_mag.begin(), reference_mag.end(), amag.begin(),
	    deltam.begin(), minus<double>());
    
  //calculate observed and abs mags in desired output bands without shift
  //calculated from SDSS bands
  k_calculate_magnitudes(coeff, redshift, zmin, zmax, band_shift,
			 filterfile, omag, amag);
  
  //apply SDSS shift to app and abs mags (same for all bands)
  vector<double>::iterator it;
  for (it=deltam.begin();it!=deltam.end();it++){
    for (i=0;i<nbands;i++){
      omag[(it-deltam.begin())*nbands+i] = omag[(it-deltam.begin())*nbands+i]+(*it);
      amag[(it-deltam.begin())*nbands+i] = amag[(it-deltam.begin())*nbands+i]+(*it);
    }
  }

  //transform(omag.begin(), omag.end(), deltam.begin(), omag.begin(), 
  //	    plus<double>());
  //transform(amag.begin(), amag.end(), deltam.begin(), amag.begin(), 
  //	    plus<double>());

}

int main(int argc, char **argv)
{
  if (argc < 2){
    cerr << "Usage: " << argv[0] << " ginfo1.dat" << endl;
    return 1;
  }

  cosmo.GetZofR(cosmo.OmegaM(),cosmo.OmegaL());
  ifstream ginfo_file(argv[1]);
  float zmin, zmax, band_shift;
  int nbands, i;
  char filterfile[FILESIZE];
  vector<ginfo> sdss_ginfo;

  if (ginfo_file.fail()) {
    cerr<<"error: cannot open "<<argv[1]<<endl;
    exit(1);
  }
  
  zmin = 0.0;
  zmax = 0.5;
  band_shift = 0.1;
  nbands = 5;
  strcpy(filterfile, "./des_filters.txt");

  copy(istream_iterator<ginfo>(ginfo_file),
	    istream_iterator<ginfo>(),
	    back_inserter(sdss_ginfo));

  int ngal = sdss_ginfo.size();
  vector<int> sed_ids(ngal);
  vector<double> redshift(ngal);
  vector<double> refmag(ngal);
  vector<double> omag(ngal*nbands);
  vector<double> amag(ngal*nbands);
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

  assign_colors(refmag, sed_ids, redshift, zmin, zmax,
		band_shift, nbands, filterfile, omag, amag);

  //write out magnitudes
  ofstream omag_file("./des_omagtest.txt");
  ofstream amag_file("./des_amagtest.txt");
  vector<double>::iterator ditr;

  for (ditr=omag.begin(); ditr!=omag.end(); ++ditr)
    {
      omag_file << *ditr;
      if ((ditr!=omag.begin()) && ((ditr-omag.begin())%nbands==0)) omag_file << "\n";
    }

  for (ditr=amag.begin(); ditr!=amag.end(); ++ditr)
    {
      amag_file << *ditr;
      if ((ditr!=amag.begin()) && ((ditr-omag.begin())%nbands==0)) amag_file << "\n";
    }

}
