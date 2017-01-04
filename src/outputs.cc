#include <vector>
#include <cassert>
#include <memory>
#include <unistd.h>     //for sleep()
#include <fstream>
#include <cstring>
#include <CCfits/CCfits>
#include "fitsio.h"
#include "outputs.h"

#ifndef UNITTESTS
int GalaxyZBin(float zRed);

void printerror( int status)
{
  /*****************************************************/
  /* Print out cfitsio error messages and exit program */
  /*****************************************************/


  if (status)
    {
      fits_report_error(stderr, status); /* print error report */

      exit( status );    /* terminate the program, returning error status */
    }
  return;
}

void PrintMinMaxDens(vector <Galaxy *> &galaxies, vector <Particle *> &particles)
{
  float min_p_dens = particles[0]->Dist8();
  float max_p_dens = particles[0]->Dist8();
  float min_g_dens = galaxies[0]->Dist8();
  float max_g_dens = galaxies[0]->Dist8();
  int min_p = 0;
  int max_p = 0;
  int min_g = 0;
  int max_g = 0;
  for(int pi=1;pi<particles.size();pi++)
    {
      if (particles[pi]->Dist8() < min_p_dens)
	{
	  min_p_dens = particles[pi]->Dist8();
	  min_p = pi;
	}
      if (particles[pi]->Dist8() > max_p_dens)
	{
	  max_p_dens = particles[pi]->Dist8();
	  max_p = pi;
	}
    }
  for(int gi=1;gi<galaxies.size();gi++)
    {
      if (galaxies[gi]->Dist8() < min_g_dens)
	{
	  min_g_dens = galaxies[gi]->Dist8();
	  min_g = gi;
	}
      if (galaxies[gi]->Dist8() > max_g_dens)
	{
	  max_g_dens = galaxies[gi]->Dist8();
	  max_g = gi;
	}
    }
  cout<<"Min/Max particle densities = "<<min_p_dens<<", "<<max_p_dens<<endl;
  cout<<"For particles "<<min_p<<", "<<max_p<<endl;
  cout<<"Min/Max galaxies densities = "<<min_g_dens<<", "<<max_g_dens<<endl;
  cout<<"For galaxies "<<min_g<<", "<<max_g<<endl;
}


void print_galaxies(vector <Galaxy *> &galaxies, vector <Particle *> &particles, vector <Halo *> &halos, vector <GalSED> galseds, vector <int> sed_ids, vector <float> nndist, vector <float> nndist_percent, string outpfn, string outdfn, string outgfn, string outghfn, string outgzfn, string outrfn)
//void print_galaxies(vector <Galaxy *> &galaxies, vector <Particle *> &particles, vector <Halo *> &halos, vector <GalSED> galseds, vector <int> sed_ids, vector <float> nndist)
{
  ofstream outpfile(outpfn.c_str());
  ofstream outdfile(outdfn.c_str());
  ofstream outgfile(outgfn.c_str());
  ofstream outghfile(outghfn.c_str());
  ofstream outgzfile(outgzfn.c_str());
  ofstream outrfile(outrfn.c_str());

  int NNotPrinted = 0;
  int NNotPrintedInVol = 0;
  for(int gi=0;gi<galaxies.size();gi++){
    Galaxy * gal = galaxies[gi];
    Particle * p = gal->P();
    assert(p);  //this better be true since you removed the other ones.
    int hid = p->Hid();
#ifdef COLORS
    int sedid = sed_ids[gi];
    if((nndist[gi]>0)&&(sedid>=0)){
      if(galseds[sed_ids[gi]].CatId() == 773)
	cout<<"Missing galaxy info = "<<gal->Mr()<<endl;
      //      outgfile<<sed_ids[gi]<<" ";
      outgfile<<galseds[sed_ids[gi]].CatId()<<" ";
      //      SEDs[sed_ids[gi]].Write(outgfile);
#else //make a dummy for the color index
    if(1){
      outgfile<<0<<" ";
#endif
      gal->Write(outgfile);
      outgzfile<<gal->zGal()<<" "<<GalaxyZBin(gal->zGal())<<" "<<gal->Central()<<endl;
      p->Write(outpfile);
      if(hid<0 || hid>=halos.size()){
	outghfile<<"0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "<<endl;
      }
      else{
	outghfile//<<in_vol<<" "
	  //<<hid<<" "
	  <<halos[hid]->Id()<<" "
	  <<halos[hid]->M()<<" "
	  <<halos[hid]->Ngal()<<" "
	  <<halos[hid]->R200()<<" "
	  <<p->Distance(halos[hid]->Position())<<" "
	  <<halos[hid]->Sig()<<" "
	  <<halos[hid]->X()<<" "
	  <<halos[hid]->Y()<<" "
	  <<halos[hid]->Z()<<" "
	  <<halos[hid]->Vx()<<" "
	  <<halos[hid]->Vy()<<" "
	  <<halos[hid]->Vz()<<" "
	  <<halos[hid]->Ra()<<" "
	  <<halos[hid]->Dec()<<" "
	  <<halos[hid]->Zred()<<endl;
	}
      outdfile<<gal->Dist8()<<" "<<gal->P()->Dist8()<<endl;
#ifdef COLORS
#ifdef COLORS_FROM_RELATIVE_DENSITY
      outrfile<<nndist[gi]<<" "<<nndist_percent[gi]<<endl;
#else
      outrfile<<nndist[gi]<<endl;
#endif
#endif
    }
    else{
      NNotPrinted++;
      if(p->Save()){
	NNotPrintedInVol++;
	cout<<"[hv] didn't print galaxy "
	    <<gi<<" "
#ifdef COLORS
	    <<nndist[gi]<<" "
	    <<sedid<<" "
#endif
	    <<p->X()<<" "
	    <<p->Y()<<" "
	    <<p->Z()<<endl;
      }
    }
  }
    cout<<"Didn't print "<<NNotPrinted<<" galaxies.  Of these "<<NNotPrintedInVol<<" were in the specified volume."<<endl;
}
#endif

void  write_bcc_catalogs(vector<Galaxy *> &galaxies, vector<Particle *> &particles,
			 vector<float> &amag, vector<float> &tmag, vector<float> &mr,
			 vector<float> &omag, vector<float> &omagerr, vector<float> &flux,
			 vector<float> &ivar, vector<float> &deltam, vector<double> &e, 
			 vector<double> &s, vector<bool> &idx, vector <Halo *> &halos,
			 vector<int> &sed_ids, vector<float> &coeffs,
			 string outgfn, string outghfn)
{
  fitsfile *fptr;
  int status, hdutype;
  long firstrow, firstelem, nrows;

  int size = galaxies.size();
  int keep = s.size();
  unsigned long rows(5);

  int tfields = 39;
  
  firstrow  = 1;
  firstelem = 1;

  char *tcolName[tfields];
  char *tcolUnit[tfields];
  char *tcolForm[tfields];

  tcolName[0] = "ID";
  tcolName[1] = "INDEX";
  tcolName[2] = "ECATID";
  tcolName[3] = "COEFFS";
  tcolName[4] = "TMAG";
  tcolName[5] = "OMAG";
  tcolName[6] = "FLUX";
  tcolName[7] = "IVAR";
  tcolName[8] = "OMAGERR";
  tcolName[9] = "AMAG";
  tcolName[10] = "RA";
  tcolName[11] = "DEC";
  tcolName[12] = "Z";
  tcolName[13] = "HALOID";
  tcolName[14] = "RHALO";
  tcolName[15] = "M200";
  tcolName[16] = "NGALS";
  tcolName[17] = "R200";
  tcolName[18] = "CENTRAL";
  tcolName[19] = "TRA";
  tcolName[20] = "TDEC";
  tcolName[21] = "EPSILON";
  tcolName[22] = "GAMMA1";
  tcolName[23] = "GAMMA2";
  tcolName[24] = "KAPPA";
  tcolName[25] = "MU";
  tcolName[26] = "MAG_R";
  tcolName[27] = "SIZE";
  tcolName[28] = "PX";
  tcolName[29] = "PY";
  tcolName[30] = "PZ";
  tcolName[31] = "VX";
  tcolName[32] = "VY";
  tcolName[33] = "VZ";
  tcolName[34] = "TE";
  tcolName[35] = "TSIZE";
  tcolName[36] = "LMAG";
  tcolName[37] = "W";
  tcolName[38] = "DELTAM";

  tcolUnit[0] = "";
  tcolUnit[1] = "";
  tcolUnit[2] = "";
  tcolUnit[3] = "";
  tcolUnit[4] = "mag";
  tcolUnit[5] = "mag";
  tcolUnit[6] = "nmgy";
  tcolUnit[7] = "nmgy^{-2}";
  tcolUnit[8] = "mag";
  tcolUnit[9] = "mag";
  tcolUnit[10] = "deg";
  tcolUnit[11] = "deg";
  tcolUnit[12] = "";
  tcolUnit[13] = "";
  tcolUnit[14] = "Mpc/h";
  tcolUnit[15] = "M_sun/h";
  tcolUnit[16] = "";
  tcolUnit[17] = "Mpc/h";
  tcolUnit[18] = "";
  tcolUnit[19] = "deg";
  tcolUnit[20] = "deg";
  tcolUnit[21] = "";
  tcolUnit[22] = "";
  tcolUnit[23] = "";
  tcolUnit[24] = "";
  tcolUnit[25] = "";
  tcolUnit[26] = "mag";
  tcolUnit[27] = "arcseconds";
  tcolUnit[28] = "Mpc/h";
  tcolUnit[29] = "Mpc/h";
  tcolUnit[30] = "Mpc/h";
  tcolUnit[31] = "km/s";
  tcolUnit[32] = "km/s";
  tcolUnit[33] = "km/s";
  tcolUnit[34] = "";
  tcolUnit[35] = "arcseconds";
  tcolUnit[36] = "mag";
  tcolUnit[37] = "";
  tcolUnit[38] = "mag";

  tcolForm[0] = "K";
  tcolForm[1] = "K";
  tcolForm[2] = "J";
  tcolForm[3] = "5E";
  tcolForm[4] = "5E";
  tcolForm[5] = "5E";
  tcolForm[6] = "5E";
  tcolForm[7] = "5E";
  tcolForm[8] = "5E";
  tcolForm[9] = "5E";
  tcolForm[10] = "E";
  tcolForm[11] = "E";
  tcolForm[12] = "E";
  tcolForm[13] = "K";
  tcolForm[14] = "E";
  tcolForm[15] = "E";
  tcolForm[16] = "J";
  tcolForm[17] = "E";
  tcolForm[18] = "J";
  tcolForm[19] = "E";
  tcolForm[20] = "E";
  tcolForm[21] = "2E";
  tcolForm[22] = "E";
  tcolForm[23] = "E";
  tcolForm[24] = "E";
  tcolForm[25] = "E";
  tcolForm[26] = "E";
  tcolForm[27] = "E";
  tcolForm[28] = "E";
  tcolForm[29] = "E";
  tcolForm[30] = "E";
  tcolForm[31] = "E";
  tcolForm[32] = "E";
  tcolForm[33] = "E";
  tcolForm[34] = "2E";
  tcolForm[35] = "E";
  tcolForm[36] = "5E";
  tcolForm[37] = "E";
  tcolForm[38] = "E";

  cout << "Extracting relevant galaxy information" << endl;
  vector<float> ra(keep), dec(keep), px(keep), py(keep), pz(keep),
    vx(keep), vy(keep), vz(keep), sdssr(keep), z(keep), id(keep),
    central(keep), haloid(keep), empty(keep,-1), cf(keep*5), am(keep*5), 
    ecatid(keep), dm(keep);

  //Go through the galaxies and get the info we want
  int count=0;
  Particle * p;
  for (int i=0; i<size; i++)
    {
#ifndef NOCUT
      if (!idx[i]) continue;
#endif
      p = galaxies[i]->P();
      int hid = p->Hid();
      ra[count] = galaxies[i]->Ra();
      dec[count] = galaxies[i]->Dec();
      id[count] = p->Gid();
      px[count] = p->X();
      py[count] = p->Y();
      pz[count] = p->Z();
      vx[count] = p->Vx();
      vy[count] = p->Vy();
      vz[count] = p->Vz();
      z[count] = p->ZredReal();
      central[count] = galaxies[i]->Central();
      sdssr[count] = mr[i];
      ecatid[count] = sed_ids[i];
      dm[count] = deltam[i];
      for (int c=0; c<5; c++)
	{
	  cf[count*5+c] = coeffs[i*5+c];
	  am[count*5+c] = amag[i*5+c];
	}
      if ((central[count]==1) & (hid>=0))
	{
	  haloid[count] = halos[hid]->Id();
	}
      else
	{
	  haloid[count] = -1;
	}
      count++;
    }

  cout << "Deleting duplicated info" << endl;
  mr.clear();
  sed_ids.clear();
  coeffs.clear();
  amag.clear();

  cout << "Number of galaxies with idx == true: " << count << endl;
  cout << "Number of galaxies with shapes: " << keep << endl;
  assert(count==keep);

  char filename[1024];
  memcpy(filename, outgfn.c_str(), outgfn.size() + 1);

  cout << "Creating file:  " << filename << endl;
  status = 0;

  fits_create_file(&fptr,filename,&status);
  if(status)
    fits_report_error(stderr,status);

  cout << "Creating TRUTH HDU" << endl;
  fits_create_tbl( fptr, BINARY_TBL, keep, tfields, tcolName, tcolForm,
		   tcolUnit, "TRUTH", &status);
  if(status)
    fits_report_error(stderr,status);      

  cout << "Getting rowsize" << endl;

  if ( fits_get_rowsize( fptr, &nrows, &status) )
    {
      cerr << "Can't get optimal number of rows." << endl;
      printerror(status);
    }

  cout << "Writing columns" << endl;

  while (firstrow<size)
    {
      fits_write_col(fptr, TLONG, 1, firstrow, firstelem, 
		     nrows, &id[0], &status);
      fits_write_col(fptr, TLONG, 2, firstrow, firstelem, 
		     nrows, &id[0], &status);
      fits_write_col(fptr, TSHORT, 3, firstrow, firstelem, 
		     nrows, &ecatid[0], &status);
      fits_write_col(fptr, TFLOAT, 4, firstrow, firstelem, 
		     nrows, &cf[0], &status);
      fits_write_col(fptr, TFLOAT, 5, firstrow, firstelem, 
		     nrows, &tmag[0], &status);
      fits_write_col(fptr, TFLOAT, 6, firstrow, firstelem, 
		     nrows, &omag[0], &status);
      fits_write_col(fptr, TFLOAT, 7, firstrow, firstelem, 
		     nrows, &flux[0], &status);
      fits_write_col(fptr, TFLOAT, 8, firstrow, firstelem, 
		     nrows, &ivar[0], &status);
      fits_write_col(fptr, TFLOAT, 9, firstrow, firstelem, 
		     nrows, &omagerr[0], &status);
      fits_write_col(fptr, TFLOAT, 10, firstrow, firstelem, 
		     nrows, &am[0], &status);
      fits_write_col(fptr, TFLOAT, 11, firstrow, firstelem, 
		     nrows, &ra[0], &status);
      fits_write_col(fptr, TFLOAT, 12, firstrow, firstelem, 
		     nrows, &dec[0], &status);
      fits_write_col(fptr, TFLOAT, 13, firstrow, firstelem, 
		     nrows, &z[0], &status);
      fits_write_col(fptr, TLONG, 14, firstrow, firstelem, 
		     nrows, &haloid[0], &status);
      fits_write_col(fptr, TFLOAT, 19, firstrow, firstelem, 
		     nrows, &central[0], &status);
      fits_write_col(fptr, TFLOAT, 22, firstrow, firstelem, 
		     nrows, &e[0], &status);
      fits_write_col(fptr, TSHORT, 27, firstrow, firstelem, 
		     nrows, &sdssr[0], &status);
      fits_write_col(fptr, TFLOAT, 28, firstrow, firstelem, 
		     nrows, &s[0], &status);
      fits_write_col(fptr, TSHORT, 29, firstrow, firstelem, 
		     nrows, &px[0], &status);
      fits_write_col(fptr, TFLOAT, 30, firstrow, firstelem, 
		     nrows, &py[0], &status);
      fits_write_col(fptr, TFLOAT, 31, firstrow, firstelem, 
		     nrows, &pz[0], &status);
      fits_write_col(fptr, TFLOAT, 32, firstrow, firstelem, 
		     nrows, &vx[0], &status);
      fits_write_col(fptr, TFLOAT, 33, firstrow, firstelem, 
		     nrows, &vy[0], &status);
      fits_write_col(fptr, TFLOAT, 34, firstrow, firstelem, 
		     nrows, &vz[0], &status);
      fits_write_col(fptr, TFLOAT, 35, firstrow, firstelem, 
		     nrows, &e[0], &status);
      fits_write_col(fptr, TFLOAT, 36, firstrow, firstelem, 
		     nrows, &s[0], &status);
      fits_write_col(fptr, TFLOAT, 40, firstrow, firstelem, 
		     nrows, &dm[0], &status);
      if (status)
	{
	  fits_report_error(stderr,status);
	  exit( status );
	}
      firstrow += nrows;
    }

}

void  write_bcc_catalogs_w_densities(
    vector<Galaxy *> &galaxies, vector<Particle *> &particles,
    vector<float> &amag, vector<float> &tmag,
    vector<float> &mr, vector<float> &omag, vector<float>
    &omagerr, vector<float> &flux, vector<float> &ivar,
    vector<float> &deltam, vector<double> &e, vector<double> &s, 
    vector<bool> &idx, vector <Halo *> &halos, vector<int> &sed_ids,
    vector<float> &coeffs, vector<float> &dist8,
    vector<float> &nndist, vector<float> &nndist_percent,
    string outgfn, string outghfn)
{

  fitsfile *fptr;
  int status, hdutype;
  long firstrow, firstelem, nrows;

  int size = galaxies.size();
  int keep = s.size();
  unsigned long rows(5);

  int tfields = 43;
  
  firstrow  = 1;
  firstelem = 1;

  char *tcolName[tfields];
  char *tcolUnit[tfields];
  char *tcolForm[tfields];

  tcolName[0] = "ID";
  tcolName[1] = "INDEX";
  tcolName[2] = "ECATID";
  tcolName[3] = "COEFFS";
  tcolName[4] = "TMAG";
  tcolName[5] = "OMAG";
  tcolName[6] = "FLUX";
  tcolName[7] = "IVAR";
  tcolName[8] = "OMAGERR";
  tcolName[9] = "AMAG";
  tcolName[10] = "RA";
  tcolName[11] = "DEC";
  tcolName[12] = "Z";
  tcolName[13] = "HALOID";
  tcolName[14] = "RHALO";
  tcolName[15] = "M200";
  tcolName[16] = "NGALS";
  tcolName[17] = "R200";
  tcolName[18] = "CENTRAL";
  tcolName[19] = "TRA";
  tcolName[20] = "TDEC";
  tcolName[21] = "EPSILON";
  tcolName[22] = "GAMMA1";
  tcolName[23] = "GAMMA2";
  tcolName[24] = "KAPPA";
  tcolName[25] = "MU";
  tcolName[26] = "MAG_R";
  tcolName[27] = "SIZE";
  tcolName[28] = "PX";
  tcolName[29] = "PY";
  tcolName[30] = "PZ";
  tcolName[31] = "VX";
  tcolName[32] = "VY";
  tcolName[33] = "VZ";
  tcolName[34] = "TE";
  tcolName[35] = "TSIZE";
  tcolName[36] = "LMAG";
  tcolName[37] = "W";
  tcolName[38] = "DIST8";
  tcolName[39] = "SIGMA5";
  tcolName[40] = "SIGMA5P";
  tcolName[41] = "PDIST8";
  tcolName[42] = "DELTAM";

  tcolUnit[0] = "";
  tcolUnit[1] = "";
  tcolUnit[2] = "";
  tcolUnit[3] = "";
  tcolUnit[4] = "mag";
  tcolUnit[5] = "mag";
  tcolUnit[6] = "nmgy";
  tcolUnit[7] = "nmgy^{-2}";
  tcolUnit[8] = "mag";
  tcolUnit[9] = "mag";
  tcolUnit[10] = "deg";
  tcolUnit[11] = "deg";
  tcolUnit[12] = "";
  tcolUnit[13] = "";
  tcolUnit[14] = "Mpc/h";
  tcolUnit[15] = "M_sun/h";
  tcolUnit[16] = "";
  tcolUnit[17] = "Mpc/h";
  tcolUnit[18] = "";
  tcolUnit[19] = "deg";
  tcolUnit[20] = "deg";
  tcolUnit[21] = "";
  tcolUnit[22] = "";
  tcolUnit[23] = "";
  tcolUnit[24] = "";
  tcolUnit[25] = "";
  tcolUnit[26] = "mag";
  tcolUnit[27] = "arcseconds";
  tcolUnit[28] = "Mpc/h";
  tcolUnit[29] = "Mpc/h";
  tcolUnit[30] = "Mpc/h";
  tcolUnit[31] = "km/s";
  tcolUnit[32] = "km/s";
  tcolUnit[33] = "km/s";
  tcolUnit[34] = "";
  tcolUnit[35] = "arcseconds";
  tcolUnit[36] = "mag";
  tcolUnit[37] = "";
  tcolUnit[38] = "Mpc/h";
  tcolUnit[39] = "Mpc/h";
  tcolUnit[40] = "";
  tcolUnit[41] = "Mpc/h";
  tcolUnit[42] = "mag";

  tcolForm[0] = "K";
  tcolForm[1] = "K";
  tcolForm[2] = "J";
  tcolForm[3] = "5E";
  tcolForm[4] = "5E";
  tcolForm[5] = "5E";
  tcolForm[6] = "5E";
  tcolForm[7] = "5E";
  tcolForm[8] = "5E";
  tcolForm[9] = "5E";
  tcolForm[10] = "E";
  tcolForm[11] = "E";
  tcolForm[12] = "E";
  tcolForm[13] = "K";
  tcolForm[14] = "E";
  tcolForm[15] = "E";
  tcolForm[16] = "J";
  tcolForm[17] = "E";
  tcolForm[18] = "J";
  tcolForm[19] = "E";
  tcolForm[20] = "E";
  tcolForm[21] = "2E";
  tcolForm[22] = "E";
  tcolForm[23] = "E";
  tcolForm[24] = "E";
  tcolForm[25] = "E";
  tcolForm[26] = "E";
  tcolForm[27] = "E";
  tcolForm[28] = "E";
  tcolForm[29] = "E";
  tcolForm[30] = "E";
  tcolForm[31] = "E";
  tcolForm[32] = "E";
  tcolForm[33] = "E";
  tcolForm[34] = "2E";
  tcolForm[35] = "E";
  tcolForm[36] = "5E";
  tcolForm[37] = "E";
  tcolForm[38] = "E";
  tcolForm[39] = "E";
  tcolForm[40] = "E";
  tcolForm[41] = "E";
  tcolForm[42] = "E";


  cout << "Extracting relevant galaxy information" << endl;
  vector<float> ra(keep), dec(keep), px(keep), py(keep), pz(keep),
    vx(keep), vy(keep), vz(keep), sdssr(keep), z(keep), id(keep),
    central(keep), empty(keep,-1), cf(keep*5),
    am(keep*5), ecatid(keep), dist8k(keep), nndistk(keep),
    nndist_percentk(keep), pdist8(keep), rhalo(keep),
    rvir(keep), mvir(keep), dm(keep);

  vector<int> hid(keep);

  //Go through the galaxies and get the info we want
  int count=0;
  Particle * p;
  for (int i=0; i<size; i++)
    {
#ifndef NOCUT
      if (!idx[i]) continue;
#endif
      p = galaxies[i]->P();
      hid[count] = p->Hid();
      rhalo[count] = p->RHalo();
      rvir[count] = p->RVir();
      mvir[count] = p->MVir();
      ra[count] = galaxies[i]->Ra();
      dec[count] = galaxies[i]->Dec();
      id[count] = p->Gid();
      px[count] = p->X();
      py[count] = p->Y();
      pz[count] = p->Z();
      vx[count] = p->Vx();
      vy[count] = p->Vy();
      vz[count] = p->Vz();
      z[count] = p->ZredReal();
      pdist8[count] = p->Dist8();
      central[count] = galaxies[i]->Central();
      sdssr[count] = mr[i];
      ecatid[count] = sed_ids[i];
      dist8k[count] = dist8[i];
      nndistk[count] = nndist[i];
      nndist_percentk[count] = nndist_percent[i];
      dm[count] = deltam[i];
      for (int c=0; c<5; c++)
	{
	  cf[count*5+c] = coeffs[i*5+c];
	  am[count*5+c] = amag[i*5+c];
	}
      if ((central[count]==1) & (hid[count]>=0))
	{
	  hid[count] = halos[hid[count]]->Id();
	}
      count++;
    }

  cout << "Deleting duplicated info" << endl;
  mr.clear();
  sed_ids.clear();
  coeffs.clear();
  amag.clear();
  dist8.clear();
  nndist.clear();
  nndist_percent.clear();
  deltam.clear();

  cout << "Number of galaxies with idx == true: " << count << endl;
  cout << "Number of galaxies with shapes: " << keep << endl;
  assert(count==keep);

  char filename[1024];
  memcpy(filename, outgfn.c_str(), outgfn.size() + 1);

  cout << "Creating file:  " << filename << endl;
  status = 0;

  fits_create_file(&fptr,filename,&status);
  if(status)
    fits_report_error(stderr,status);

  cout << "Creating TRUTH HDU" << endl;
  fits_create_tbl( fptr, BINARY_TBL, keep, tfields, tcolName, tcolForm,
		   tcolUnit, "TRUTH", &status);
  if(status)
    fits_report_error(stderr,status);      

  cout << "Getting rowsize" << endl;

  if ( fits_get_rowsize( fptr, &nrows, &status) )
    {
      cerr << "Can't get optimal number of rows." << endl;
      printerror(status);
    }

  cout << "Writing columns" << endl;

  while (firstrow<size)
    {
      fits_write_col(fptr, TLONG, 1, firstrow, firstelem, 
		     nrows, &id[0], &status);
      fits_write_col(fptr, TLONG, 2, firstrow, firstelem, 
		     nrows, &id[0], &status);
      fits_write_col(fptr, TSHORT, 3, firstrow, firstelem, 
		     nrows, &ecatid[0], &status);
      fits_write_col(fptr, TFLOAT, 4, firstrow, firstelem, 
		     nrows, &cf[0], &status);
      fits_write_col(fptr, TFLOAT, 5, firstrow, firstelem, 
		     nrows, &tmag[0], &status);
      fits_write_col(fptr, TFLOAT, 6, firstrow, firstelem, 
		     nrows, &omag[0], &status);
      fits_write_col(fptr, TFLOAT, 7, firstrow, firstelem, 
		     nrows, &flux[0], &status);
      fits_write_col(fptr, TFLOAT, 8, firstrow, firstelem, 
		     nrows, &ivar[0], &status);
      fits_write_col(fptr, TFLOAT, 9, firstrow, firstelem, 
		     nrows, &omagerr[0], &status);
      fits_write_col(fptr, TFLOAT, 10, firstrow, firstelem, 
		     nrows, &am[0], &status);
      fits_write_col(fptr, TFLOAT, 11, firstrow, firstelem, 
		     nrows, &ra[0], &status);
      fits_write_col(fptr, TFLOAT, 12, firstrow, firstelem, 
		     nrows, &dec[0], &status);
      fits_write_col(fptr, TFLOAT, 13, firstrow, firstelem, 
		     nrows, &z[0], &status);
      fits_write_col(fptr, TLONG, 14, firstrow, firstelem, 
		     nrows, &hid[0], &status);
      fits_write_col(fptr, TFLOAT, 19, firstrow, firstelem, 
		     nrows, &central[0], &status);
      fits_write_col(fptr, TFLOAT, 22, firstrow, firstelem, 
		     nrows, &e[0], &status);
      fits_write_col(fptr, TSHORT, 27, firstrow, firstelem, 
		     nrows, &sdssr[0], &status);
      fits_write_col(fptr, TFLOAT, 28, firstrow, firstelem, 
		     nrows, &s[0], &status);
      fits_write_col(fptr, TSHORT, 29, firstrow, firstelem, 
		     nrows, &px[0], &status);
      fits_write_col(fptr, TFLOAT, 30, firstrow, firstelem, 
		     nrows, &py[0], &status);
      fits_write_col(fptr, TFLOAT, 31, firstrow, firstelem, 
		     nrows, &pz[0], &status);
      fits_write_col(fptr, TFLOAT, 32, firstrow, firstelem, 
		     nrows, &vx[0], &status);
      fits_write_col(fptr, TFLOAT, 33, firstrow, firstelem, 
		     nrows, &vy[0], &status);
      fits_write_col(fptr, TFLOAT, 34, firstrow, firstelem, 
		     nrows, &vz[0], &status);
      fits_write_col(fptr, TFLOAT, 35, firstrow, firstelem, 
		     nrows, &e[0], &status);
      fits_write_col(fptr, TFLOAT, 36, firstrow, firstelem, 
		     nrows, &s[0], &status);
      fits_write_col(fptr, TFLOAT, 39, firstrow, firstelem, 
		     nrows, &dist8k[0], &status);
      fits_write_col(fptr, TFLOAT, 40, firstrow, firstelem, 
		     nrows, &nndistk[0], &status);
      fits_write_col(fptr, TFLOAT, 41, firstrow, firstelem, 
		     nrows, &nndist_percentk[0], &status);
      fits_write_col(fptr, TFLOAT, 42, firstrow, firstelem, 
		     nrows, &pdist8[0], &status);
      fits_write_col(fptr, TFLOAT, 43, firstrow, firstelem, 
		     nrows, &dm[0], &status);
      if (status)
	{
	  fits_report_error(stderr,status);
	  exit( status );
	}
      firstrow += nrows;
    }
}
