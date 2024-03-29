#include "kcorrect_utils.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <iomanip>
#include <algorithm>
#define FILESIZE 2000
#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}

using namespace std;

static int nz=1500;
static int nk,nv,nl,maxn;
static float *lambda=NULL;
static float *vmatrix=NULL;
static float *rmatrix=NULL;
static float *zvals=NULL;
static float *filter_lambda=NULL;
static float *filter_pass=NULL;
static int *filter_n=NULL;

istream & operator>>(istream & is, magtuple & in)
{
  is >> in.bands[0] >> in.bands[1] >> in.bands[2] >>
    in.bands[3] >> in.bands[4];

  return is;
}

void reconstruct_maggies(float *coeff, float *redshift, int ngal, float zmin, 
			 float zmax, float band_shift, char filterfile[], float *maggies)
{
  cout<<"Reconstructing maggies"<<endl;

  int i,ndim,*sizes=NULL;
  char vfile[FILESIZE],lfile[FILESIZE],path[FILESIZE];
  char vmatrixfile[FILESIZE],lambdafile[FILESIZE];

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
  cout<<"Loading filters from "<<filterfile<<endl;
  k_load_filters(&filter_n,&filter_lambda,&filter_pass,&maxn,&nk,filterfile);
  cout<<"Number of filters in "<<filterfile<<": "<<nk<<endl;

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

void match_coeff(vector<int> &sed_ids, float* coeffs)
{
  static int ntemps=5;
  int n;
  //Need to replace with variable passed in stringparameters
  string filename = "/nfs/slac/g/ki/ki23/des/jderose/l-addgals/training/cooper/combined_dr6_cooper.dat";
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

  for (n=0;n<ntemps;n++){
    cout<<setprecision(numeric_limits<float>::digits10+1)
	<<sdsscoeffs[0].bands[n]<<", ";
  }
  cout<<endl;
  
  for (it = sed_ids.begin(); it!=sed_ids.end(); it++)
    {
      for (n=0; n<ntemps; ++n) {
	coeffs[(it-sed_ids.begin())*ntemps+n] = sdsscoeffs[*it-1].bands[n];
      }
    }
}

void k_calculate_magnitudes(vector<float> &coeff, vector<float> &redshift,
			    float zmin, float zmax, float band_shift,
			    int nband, char filterfile[], vector<float> &omag,
			    vector<float> &amag)
{
  cout<<"Calculating magnitudes"<<endl;
  int i;
  int ngal = redshift.size();
  float zeropoint = 22.5;
  vector<float> kcorrection(omag.size());
  vector<float> rf_z(ngal);
  fill(rf_z.begin(), rf_z.end(), band_shift);

  //observer frame maggies
  reconstruct_maggies(&coeff[0], &redshift[0], ngal, zmin,
                      zmax, 0.0, filterfile, &omag[0]);

  //rest frame maggies
  reconstruct_maggies(&coeff[0], &rf_z[0], ngal, zmin,
                      zmax, 0.0, filterfile, &amag[0]);

  cout<<"Calculating kcorrections"<<endl;
  //calculate k-correction in mags
  transform(amag.begin(), amag.end(), omag.begin(), 
	    kcorrection.begin(), magop<float>());

  cout<<"First few kcorrections: "<<endl;
  for (i=0; i<nband; i++){
    cout<<kcorrection[i]<<", ";
  }
  cout<<endl;

  cout<<"Calculating apparent magnitudes"<<endl;
  //calculate apparent magnitudes
  transform(omag.begin(), omag.end(), omag.begin(),
	    bind2nd(appmagnitude<float>(),zeropoint));

  cout<<"First few apparent mags: "<<endl;
  for (i=0; i<nband; i++){
    cout<<omag[i]<<", ";
  }
  cout<<endl;

  cout<<"Calculating absolute magnitudes"<<endl;
  //calculate absolute magnitudes without kcorrection
  
  vector<float>::iterator it;
  absmagnitude<float> absmagn;
  for (it=redshift.begin(); it!=redshift.end(); it++){
    for (i=0; i<nband; i++){
      amag[nband*(it-redshift.begin())+i] = absmagn(omag[nband*(it-redshift.begin())+i],*it);
    }
  }

  cout<<"Subtracting kcorrection"<<endl;
  //add k-correction to absolute magnitudes
  transform(amag.begin(), amag.end(), kcorrection.begin(),
	    amag.begin(), minus<float>());

  cout<<"First few absolute mags: "<<endl;
  for (i=0; i<nband; i++){
    cout<<amag[i]<<", ";
  }
  cout<<endl;

  cout<<"Done calculating magnitudes"<<endl;
}

void assign_colors(vector<float> &reference_mag, vector<float> &coeff, 
		   vector<float> &redshift, float zmin, float zmax,
		   float band_shift, int nbands, char filterfile[], 
		   vector<float> &omag, vector<float> &amag){
  int i;
  int ngal = reference_mag.size();
  int ntemp = 5;
  float sdss_bandshift = 0.1;
  float ab_corr = 0.012;
  char sdss_filterfile[FILESIZE];
  vector<float> deltam(ngal);
  strcpy(sdss_filterfile, "/nfs/slac/g/ki/ki23/des/jderose/l-addgals/src/sdss_filter.txt");
  
  cout << "First few reference magnitudes: " << endl;
  for (i=0;i<5;i++) 
    {
      cout << reference_mag[i] << ", ";
    }
  cout << endl;

  cout<<"Coeffs: "<<endl;
  for (i=0;i<25;i++){
    cout<<coeff[i]<<", ";
    if ((i+1)%5==0) cout<<endl;
  }
  cout<<endl;

  //calculate SDSS r-band abs mags to determine magnitude shifts to apply
  //to kcorrect outputs for other bands
  k_calculate_magnitudes(coeff, redshift, zmin, zmax, sdss_bandshift,
			 1, sdss_filterfile, omag, amag);

  //apply ab magnitude correction
  transform(amag.begin(), amag.end(), amag.begin(),
  	    bind2nd(minus<float>(), ab_corr));

  //determine shift needed to get reference_mag from sdss_amag
  transform(reference_mag.begin(), reference_mag.end(), amag.begin(),
	    deltam.begin(), minus<float>());

  cout<<"10-100 deltam"<<endl;
  for (i=10;i<100;i++){
    cout<<deltam[i]<<", ";
    if ((i-10+1)%5==0) cout<<endl;
  }
  cout<<endl;
  
  //calculate observed and abs mags in desired output bands without shift
  //calculated from SDSS bands
  k_calculate_magnitudes(coeff, redshift, zmin, zmax, band_shift,
			 nbands, filterfile, omag, amag);
  
  //apply SDSS shift to app and abs mags (same for all bands)
  vector<float>::iterator it;
  for (it=deltam.begin();it!=deltam.end();it++){
    for (i=0;i<nbands;i++){
      omag[(it-deltam.begin())*nbands+i] = omag[(it-deltam.begin())*nbands+i]+(*it);
      amag[(it-deltam.begin())*nbands+i] = amag[(it-deltam.begin())*nbands+i]+(*it);
    }
  }

#ifdef COLORTEST
  ofstream amag_file("./des_amagtest.txt");
  vector<float>::iterator ditr;

  for (ditr=amag.begin(); ditr!=amag.end(); ++ditr)
    {
      amag_file << *ditr;
      if ((ditr-amag.begin()+1)%nbands==0)
	{
	  amag_file << "\n";
	} else amag_file << " ";
    }
#endif

}

#ifdef UNITTESTS

int main(int argc, char **argv)
{
  if (argc < 2){
    cerr << "Usage: " << argv[0] << " ginfo1.dat" << endl;
    return 1;
  }

  ifstream ginfo_file(argv[1]);
  float zmin, zmax, band_shift;
  int nbands, ntemp, i;
  char filterfile[FILESIZE];
  float ZREDMIN = 0.00000;
  float ZREDMAX = 0.174762;

  vector<ginfo> sdss_ginfo;

  if (ginfo_file.fail()) {
    cerr<<"error: cannot open "<<argv[1]<<endl;
    exit(1);
  }
  zmin = 0.0;
  zmax = 2.5;
  band_shift = 0.1;
  nbands = 5;
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
  ofstream cf_file("./des_cftest.txt");
  vector<float>::iterator ditr;

  for (ditr=omag.begin(); ditr!=omag.end(); ++ditr)
    {
      omag_file << *ditr;
      if ((ditr-omag.begin()+1)%nbands==0)
	{
	  omag_file << "\n";
	} else omag_file << " ";
    }

  for (ditr=coeff.begin(); ditr!=coeff.end(); ++ditr)
    {
      cf_file << *ditr;
      if ((ditr-coeff.begin()+1)%nbands==0)
	{
	  cf_file << "\n";
	} else cf_file << " ";
    }

  for (ditr=amag.begin(); ditr!=amag.end(); ++ditr)
    {
      amag_file << *ditr;
      if ((ditr-amag.begin()+1)%nbands==0)
	{
	  amag_file << "\n";
	} else amag_file << " ";
    }
}
#endif
