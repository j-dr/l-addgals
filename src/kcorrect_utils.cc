#include "kcorrect_utils.h"
#include "hv.h"
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <stdlib.h>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <iomanip>
#include <algorithm>
#include <CCfits/CCfits>
#include <valarray>
#include "fitsio.h"
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

#ifdef MAKEMAGS
string colortrdir;
#endif
#ifdef UNITTESTS
string colortrdir = "/nfs/slac/g/ki/ki23/des/jderose/SkyFactory-config/Addgals";
#endif

void reconstruct_maggies(float *coeff, float *redshift, int ngal, float zmin_this,
			 float zmax_this, float band_shift, char filterfile[], float *maggies)
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
    zvals[i]=zmin_this+(zmax_this-zmin_this)*((float)i+0.5)/(float)nz;
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
  string filename = colortrdir + "/combined_dr6_cooper.dat";
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
			    float zmin_this, float zmax_this, float band_shift,
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
  reconstruct_maggies(&coeff[0], &redshift[0], ngal, zmin_this,
		      zmax_this, 0.0, filterfile, &omag[0]);

  //rest frame maggies
  reconstruct_maggies(&coeff[0], &rf_z[0], ngal, zmin_this,
		      zmax_this, 0.0, filterfile, &amag[0]);

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
		   vector<float> &redshift, float zmin_this, float zmax_this,
		   float band_shift, int nbands, char filterfile[],
		   vector<float> &omag, vector<float> &amag, vector<float> &deltam,
		   vector<float> abcorr, bool refflag){

  int i;
  int ngal = reference_mag.size();
  int ntemp = 5;
  float sdss_bandshift = 0.1;
  float sdssab_corr = 0.012;
  char sdss_filterfile[FILESIZE];
  strcpy(sdss_filterfile, (colortrdir+"/sdss_filter.txt").c_str());

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
  k_calculate_magnitudes(coeff, redshift, zmin_this, zmax_this, sdss_bandshift,
			 1, sdss_filterfile, omag, amag);

  //determine shift needed to get reference_mag from sdss_amag
  if (refflag)
    {
      transform(amag.begin(), amag.end(), amag.begin(),
		bind2nd(minus<float>(), sdssab_corr));
      transform(reference_mag.begin(), reference_mag.end(), amag.begin(),
		deltam.begin(), minus<float>());
    }
  else
    {
      transform(omag.begin(), omag.end(), omag.begin(),
		bind2nd(minus<float>(), sdssab_corr));

      transform(reference_mag.begin(), reference_mag.end(), omag.begin(),
		deltam.begin(), minus<float>());
    }

  cout<<"10-100 deltam"<<endl;
  for (i=10;i<100;i++){
    cout<<deltam[i]<<", ";
    if ((i-10+1)%5==0) cout<<endl;
  }
  cout<<endl;

  //calculate observed and abs mags in desired output bands without shift
  //calculated from SDSS bands
  k_calculate_magnitudes(coeff, redshift, zmin_this, zmax_this, band_shift,
			 nbands, filterfile, omag, amag);

  //apply SDSS shift to app and abs mags (same for all bands)
  vector<float>::iterator it;
  for (it=deltam.begin();it!=deltam.end();it++){
    for (i=0;i<nbands;i++){
      omag[(it-deltam.begin())*nbands+i] = omag[(it-deltam.begin())*nbands+i]+(*it)+abcorr[i];
      amag[(it-deltam.begin())*nbands+i] = amag[(it-deltam.begin())*nbands+i]+(*it)+abcorr[i];
    }
    for (i=0;i<ntemp;i++)
      coeff[(it-deltam.begin())*ntemp+i] *= pow(10.0, *it/-2.5);
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

long readNRowsFits(std::string filename)
{
  using namespace CCfits;

  std::vector<std::string> hdukeys(2,"");
  hdukeys[0] = "MAG_R";
  hdukeys[1] = "COEFFS";

  //Create fits object
  std::auto_ptr<FITS> pInfile(new FITS(filename, Read, 1, false, hdukeys));

  ExtHDU& table = pInfile->extension(1);

  std::cout << "NAXIS2: " << table.rows() << std::endl;

  return table.rows();
}

//void readSEDInfoFITS(std::string filename, std::vector< valarray<float> > &coeffs, std::vector<float> &sdss_mag_r)
void readSEDInfoFITS(std::string filename, std::vector<float> &coeffs, std::vector<float> &sdss_mag_r, std::vector<float> &redshift, std::vector<int> &id)
{
  using namespace CCfits;

  const long nrows = sdss_mag_r.size();
  long i;
  std::vector< std::valarray<float> > temp(nrows);

  //Create fits object
  std::auto_ptr<FITS> pInfile(new FITS(filename, Read, 1));
  ExtHDU& table = pInfile->extension(1);

  //read sdss mag, kcorrect coeffs, and redshifts
  table.column("MAG_R").read(sdss_mag_r,1,nrows);
  table.column("Z").read(redshift,1,nrows);
  table.column("ECATID").read(id,1,nrows);

  vector<float>::iterator itr=coeffs.begin();
  vector< valarray<float> >::iterator titr;
  table.column("COEFFS").readArrays(temp, 1, nrows);
  for (titr=temp.begin();titr<temp.end();titr++)
    {
      for (i=0;i<5;i++)
	{
	  coeffs[(titr-temp.begin())*5+i] = (*titr)[i];
	}
    }
}

void write_colors(std::vector<float> amag, std::vector<float> omag, std::vector<float> coeff_norm,
		  const int nbands, std::string outname, std::string surveyname)
{

  using namespace CCfits;
  std::auto_ptr<FITS> tFits;
  int size = amag.size()/nbands;
  unsigned long rows(nbands);

  //Write truth file
  cout << "Opening fits file" << endl;
  try{
    tFits.reset(new FITS(outname,Write));
  }
  catch (CCfits::FITS::CantOpen){
    cerr << "Can't open " << outname << endl;
  }

  vector<string> tcolName(6,"");
  vector<string> tcolUnit(6,"");
  vector<string> tcolForm(6,"");

  tcolName[0] = "AMAG";
  tcolName[1] = "TMAG";
  tcolName[2] = "OMAG";
  tcolName[3] = "OMAGERR";
  tcolName[4] = "FLUX";
  tcolName[5] = "IVAR";

  tcolUnit[0] = "mag";
  tcolUnit[1] = "mag";
  tcolUnit[2] = "mag";
  tcolUnit[3] = "mag";
  tcolUnit[4] = "nmgy";
  tcolUnit[5] = "nmgy^{-2}";


  std::string cf;
  std::stringstream strstream;
  strstream << nk << "E";
  strstream >> cf;

  tcolForm[0] = cf;
  tcolForm[1] = cf;
  tcolForm[2] = cf;
  tcolForm[3] = cf;
  tcolForm[4] = cf;
  tcolForm[5] = cf;

  Table* newTable;

  try{
    newTable = tFits->addTable(surveyname,size,tcolName,tcolForm,tcolUnit);
  }
  catch(...){
    printf("Could not create table\n");
    exit(1);
  }

  try{
    newTable->column(tcolName[0]).write(amag,size,1);
    newTable->column(tcolName[1]).write(omag,size,1);
  }
  catch(FitsException &except){
    printf("Caught Save Error: Column Write -- ");
    printf("%s\n",except.message().c_str());
    exit(1);
  }
  catch(...){
    printf("Caught Save Error: Column Write\n");
    exit(1);
  }


}



void write_colors_cfitsio(std::vector<float> amag, std::vector<float> omag, std::vector<float> coeff_norm,
		  const int nbands, std::string outname, std::string surveyname)
{
  
  fitsfile *fptr;
  int status, hdutype;
  long firstrow, firstelem;

  int tfields = 6;

  char filename[1024];
  //char filename[] = "test.fits";
  memcpy(filename, outname.c_str(), outname.size() + 1);

  char extname[]  = "SURVEYMAGS";

  std::string cf;
  std::stringstream strstream;
  strstream << nbands << "E";
  strstream >> cf;

  char *cff =  (char *)malloc(cf.size() + 1);
  memcpy(cff, cf.c_str(), cf.size() + 1);

  char *ttype[]    = {"AMAG", "TMAG", "OMAG", "OMAGERR", "FLUX", "IVAR"};
  char *tcolunit[] = {"mag", "mag", "mag", "mag", "nmgy", "nmgy^{-2}"};
  char *tform[]   = {cff, cff, cff, cff, cff, cff};

  int size = amag.size()/nbands;
  long nrows;

  //Write truth file

  cout << "Creating file:  " << filename << endl;
  status = 0;

  fits_create_file(&fptr,filename,&status);
  if(status)
    fits_report_error(stderr,status);

  fits_create_tbl( fptr, BINARY_TBL, size, tfields, ttype, tform,
		   tcolunit, extname, &status);
  if(status)
    fits_report_error(stderr,status);      

  firstrow  = 1;
  firstelem = 1;

  cout << "Getting rowsize" << endl;

  if ( fits_get_rowsize( fptr, &nrows, &status) )
    {
      cerr << "Can't get optimal number of rows." << endl;
      fits_report_error(stderr, status);
    }

  cout << "Writing columns" << endl;

  while (firstrow<size)
    {
      fits_write_col(fptr, TFLOAT, 1, firstrow, firstelem, 
		     nrows, &amag[0], &status);
      fits_write_col(fptr, TFLOAT, 2, firstrow, firstelem, 
		     nrows, &omag[0], &status);
      firstrow += nrows;
    }

  cout << "Closing file" << endl;

  if (fits_close_file(fptr, &status) )
    {
      cerr << "Can't close file " << outname << endl;
      cerr << "status: close file " << status << endl;
    }

}


#ifdef MAKEMAGS

enum string_code
  {
  TWOMASS,
  CFHTLS,
  CANDELS,
  DEEP2,
  SDSS,
  DECAM,
  Euclid,
  FLAMEX,
  HSC,
  IRAC,
  Johnson,
  LSST,
  VISTA,
  WFIRST,
  WISE,
  DESY3
};

string_code hashit (std::string const& inString)
{
  if (inString == "TWOMASS") return TWOMASS;
  if (inString == "CFHTLS") return CFHTLS;
  if (inString == "CANDELS") return CANDELS;
  if (inString == "DEEP2") return DEEP2;
  if (inString == "SDSS") return SDSS;
  if (inString == "DECAM") return DECAM;
  if (inString == "DESY3") return DESY3;
  if (inString == "Euclid") return Euclid;
  if (inString == "FLAMEX") return FLAMEX;
  if (inString == "HSC") return HSC;
  if (inString == "IRAC") return IRAC;
  if (inString == "Johnson") return Johnson;
  if (inString == "LSST") return LSST;
  if (inString == "VISTA") return VISTA;
  if (inString == "WFIRST") return WFIRST;
  if (inString == "WISE") return WISE;
}

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
  std::vector<std::string> elems;
  split(s, delim, elems);
  return elems;
}

int main(int argc, char *argv[])
{
  long nrows;
  int i,j;
  int ntemp = 5;
  float band_shift=0.0;
  float t1, t2;

  std::string filename(argv[1]);
  std::string outdir(argv[2]);
  std::string outbase(argv[3]);
  colortrdir = argv[4];

  cout << "Input file name  : " << filename << endl;
  cout << "Output directory : " << outdir << endl;
  cout << "Outbase          : " << outbase << endl;
  cout << "Colortrdir       : " << colortrdir << endl;

  std::string survey;
  std::string outname;
  char filterfile[FILESIZE];

  std::vector<std::string> psplt = split(filename, '.');

  //Read in galaxy SED information
  nrows = readNRowsFits(filename);

  std::vector<float> coeffs(nrows*ntemp);
  std::vector<float> sdss_mag_r(nrows);
  std::vector<float> redshift(nrows);
  std::vector<int> id(nrows);

  readSEDInfoFITS(filename, coeffs, sdss_mag_r, redshift, id);

  match_coeff(id, &coeffs[0]);
  //Loop over surveys
  for (i=5;i<argc;i++)
    {
      survey = argv[i];
      cout << "Survey      : " << survey << endl;
      //Load filter information
      float aabcorr[5] = {-0.036, 0.012, 0.010, 0.028, 0.040};
      vector<float> abcorr;

      switch(hashit(survey))
	{
	case TWOMASS:
	  strcpy(filterfile, (colortrdir+"/twomass_filters.txt").c_str());
	  break;
	case CFHTLS:
	  strcpy(filterfile, (colortrdir+"/cfhtls_filters.txt").c_str());
	  break;
	case CANDELS:
	  strcpy(filterfile, (colortrdir+"/candels_filters.txt").c_str());
	  break;
	case DEEP2:
	  strcpy(filterfile, (colortrdir+"/deep2_filters.txt").c_str());
	  break;
	case SDSS:
	  strcpy(filterfile, (colortrdir+"/sdss_filters.txt").c_str());
	  band_shift = 0.1;
	  break;
	case DECAM:
	  strcpy(filterfile, (colortrdir+"/decam_filters.txt").c_str());
	  band_shift = 0.1;
	  break;
        case DESY3:
	  strcpy(filterfile, (colortrdir+"des_y3_filters.txt").c_str());
	  band_shift = 0.1;
	  break;
	case Euclid:
	  strcpy(filterfile, (colortrdir+"/euclid_filters.txt").c_str());
	  break;
	case FLAMEX:
	  strcpy(filterfile, (colortrdir+"/flamex_filters.txt").c_str());
	  break;
	case HSC:
	  strcpy(filterfile, (colortrdir+"/hsc_filters.txt").c_str());
	  break;
	case IRAC:
	  strcpy(filterfile, (colortrdir+"/irac_filters.txt").c_str());
	  break;
	case Johnson:
	  strcpy(filterfile, (colortrdir+"/johnson_filters.txt").c_str());
	  break;
	case LSST:
	  strcpy(filterfile, (colortrdir+"/lsst_filters.txt").c_str());
	  break;
	case VISTA:
	  strcpy(filterfile, (colortrdir+"/vista_filters.txt").c_str());
	  band_shift = 0.1;
	  break;
	case WFIRST:
	  strcpy(filterfile, (colortrdir+"/wfirst_filters.txt").c_str());
	  break;
	case WISE:
	  strcpy(filterfile, (colortrdir+"/wise_filters.txt").c_str());
	  break;
	}

      cout<<"Loading filters from "<<filterfile<<endl;
      k_load_filters(&filter_n,&filter_lambda,&filter_pass,&maxn,&nk,filterfile);
      cout<<"Number of filters in "<<filterfile<<": "<<nk<<endl;

      std::vector<float> omag(nrows*nk);
      std::vector<float> amag(nrows*nk);
      std::vector<float> coeff_norm(nrows);

      //generate survey magnitudes
      if (hashit(survey)==SDSS)
	{
	  for (j=0;j<nk;j++)
	    {
	      abcorr.push_back(aabcorr[j]);
	    }

	  assign_colors(sdss_mag_r, coeffs, redshift, 0.0, 2.5, band_shift,
			nk, filterfile, omag, amag, coeff_norm, abcorr);
	}
      else
	{
	  for (j=0;j<nk;j++)
	    {
	      abcorr.push_back(0.0);
	    }

      	  assign_colors(sdss_mag_r, coeffs, redshift, 0.0, 2.5, band_shift,
			nk, filterfile, omag, amag, coeff_norm, abcorr);
	}
      
      outname = outdir + "/" + outbase + "_" + survey + "." + *(psplt.end()-2) + ".fits";
      cout << "Writing to " << outname << endl;

      t1 = clock();      
      write_colors_cfitsio(amag, omag, coeff_norm, nk, outname, survey);
      t2 = clock();
      cout << "cfitsio wrote outputs in " << (t2-t1)/CLOCKS_PER_SEC << "seconds" << endl;      

      //      t1 = clock();      
      //      write_colors(amag, omag, coeff_norm, nk, outname, survey);
      //      t2 = clock();
      //      cout << "ccfits wrote outputs in " << (t2-t1)/CLOCKS_PER_SEC << "seconds" << endl;      
    }
}

#endif

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
  vector<float> coeff_norm(ngal);
  vector<float> omag(ngal*nbands);
  vector<float> amag(ngal*nbands);
  vector<float> coeff(ngal*ntemp);
  vector<float> abcorr(5,0.0);
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
		band_shift, nbands, filterfile, omag,
		amag,coeff_norm,abcorr,false);

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
