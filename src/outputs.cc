#include <vector>
#include <cassert>
#include <memory>
#include <unistd.h>     //for sleep()
#include <fstream>
#include <CCfits/CCfits>
#include "outputs.h"

#ifndef UNITTESTS
int GalaxyZBin(float zRed);

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
			 vector<float> &ivar, vector<double> &e, vector<double> &s,
			 vector<bool> &idx, vector <Halo *> &halos,
			 vector<int> &sed_ids, vector<float> &coeffs,
			 string outgfn, string outghfn)
{
  using namespace CCfits;
  std::auto_ptr<FITS> tFits;
  int size = galaxies.size();
  int keep = s.size();
  unsigned long rows(5);

  //Write truth file
  cout << "Opening fits file" << endl;
  try{
    tFits.reset(new FITS(outgfn,Write));
  }
  catch (CCfits::FITS::CantOpen){
    cerr << "Can't open " << outgfn << endl;
  }

  vector<string> tcolName(36,"");
  vector<string> tcolUnit(36,"");
  vector<string> tcolForm(36,"");

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


  cout << "Extracting relevant galaxy information" << endl;
  vector<float> ra(keep), dec(keep), px(keep), py(keep), pz(keep),
    vx(keep), vy(keep), vz(keep), sdssr(keep), z(keep), id(keep),
    central(keep), haloid(keep), empty(keep,-1), cf(keep*5), am(keep*5), ecatid(keep);

  //Go through the galaxies and get the info we want
  int count=0;
  Particle * p;
  for (int i=0; i<size; i++)
    {
      if (!idx[i]) continue;
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
      z[count] = p->Zred();
      central[count] = galaxies[i]->Central();
      sdssr[count] = mr[i];
      ecatid[count] = sed_ids[i];
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

  cout << "Number of galaxies with idx == true: " << count << endl;
  cout << "Number of galaxies with shapes: " << keep << endl;
  assert(count==keep);

  cout << "Creating TRUTH HDU" << endl;
  static string ttablename("TRUTH");
  Table* newTable;

  try{
    newTable = tFits->addTable(ttablename,keep,tcolName,tcolForm,tcolUnit);
  }
  catch(...){
    printf("Could not create table\n");
    exit(1);
  }

  cout << "Writing columns" << endl;
  try{
    cout << "writing id 1" << endl;
    newTable->column(tcolName[0]).write(id,1);
    cout << "writing id 2" << endl;
    newTable->column(tcolName[1]).write(id,1);
    cout << "writing sed_ids" << endl;
    newTable->column(tcolName[2]).write(ecatid,1);
    cout << "writing coeffs" << endl;
    newTable->column(tcolName[3]).write(cf,count,1);
    cout << "writing tmag" << endl;
    newTable->column(tcolName[4]).write(tmag,count,1);
    cout << "writing omag" << endl;
    newTable->column(tcolName[5]).write(omag,count,1);
    cout << "writing omag" << endl;
    newTable->column(tcolName[6]).write(flux,count,1);
    cout << "writing flux" << endl;
    newTable->column(tcolName[7]).write(ivar,count,1);
    cout << "writing ivar" << endl;
    newTable->column(tcolName[8]).write(omagerr,count,1);
    cout << "writing coeffs" << endl;
    newTable->column(tcolName[9]).write(am,count,1);
    cout << "writing ra" << endl;
    newTable->column(tcolName[10]).write(ra,1);
    cout << "writing dec" << endl;
    newTable->column(tcolName[11]).write(dec,1);
    cout << "writing z" << endl;
    newTable->column(tcolName[12]).write(z,1);
    cout << "writing haloid" << endl;
    newTable->column(tcolName[13]).write(haloid,1);
    cout << "writing 18" << endl;
    newTable->column(tcolName[18]).write(central,1);
    cout << "writing 21" << endl;
    newTable->column(tcolName[21]).write(e,count,1);
    cout << "writing 26" << endl;
    newTable->column(tcolName[26]).write(sdssr,1);
    cout << "writing 27" << endl;
    newTable->column(tcolName[27]).write(s,1);
    cout << "writing 28" << endl;
    newTable->column(tcolName[28]).write(px,1);
    cout << "writing 29" << endl;
    newTable->column(tcolName[29]).write(py,1);
    cout << "writing 30" << endl;
    newTable->column(tcolName[30]).write(pz,1);
    cout << "writing 31" << endl;
    newTable->column(tcolName[31]).write(vx,1);
    cout << "writing 32" << endl;
    newTable->column(tcolName[32]).write(vy,1);
    cout << "writing 33" << endl;
    newTable->column(tcolName[33]).write(vz,1);
    cout << "writing 34" << endl;
    newTable->column(tcolName[34]).write(e,count,1);
    cout << "writing 35" << endl;
    newTable->column(tcolName[35]).write(s,1);
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

void  write_bcc_catalogs_w_densities(vector<Galaxy *> &galaxies, vector<Particle *> &particles,
			 vector<float> &amag, vector<float> &tmag, vector<float> &mr,
			 vector<float> &omag, vector<float> &omagerr, vector<float> &flux,
			 vector<float> &ivar, vector<double> &e, vector<double> &s,
			 vector<bool> &idx, vector <Halo *> &halos,
			 vector<int> &sed_ids, vector<float> &coeffs,
			 vector<float> &dist8, vector<float> &nndist,
             vector<float> &nndist_percent,
             vector<float> &string outgfn, string outghfn)
{
  using namespace CCfits;
  std::auto_ptr<FITS> tFits;
  int size = galaxies.size();
  int keep = s.size();
  unsigned long rows(5);

  //Write truth file
  cout << "Opening fits file" << endl;
  try{
    tFits.reset(new FITS(outgfn,Write));
  }
  catch (CCfits::FITS::CantOpen){
    cerr << "Can't open " << outgfn << endl;
  }

  vector<string> tcolName(36,"");
  vector<string> tcolUnit(36,"");
  vector<string> tcolForm(36,"");

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
  tcolName[38] = "DIST8"
  tcolName[39] = "SIGMA5"
  tcolName[40] = "SIGMA5P"

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
  tcolUnit[38] = "Mpc/h"
  tcolUnit[39] = "Mpc/h"
  tcolUnit[40] = ""

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


  cout << "Extracting relevant galaxy information" << endl;
  vector<float> ra(keep), dec(keep), px(keep), py(keep), pz(keep),
    vx(keep), vy(keep), vz(keep), sdssr(keep), z(keep), id(keep),
    central(keep), haloid(keep), empty(keep,-1), cf(keep*5),
    am(keep*5), ecatid(keep), dist8k(keep), nndistk(keep),
    nndist_percentk(keep);

  //Go through the galaxies and get the info we want
  int count=0;
  Particle * p;
  for (int i=0; i<size; i++)
    {
      if (!idx[i]) continue;
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
      z[count] = p->Zred();
      central[count] = galaxies[i]->Central();
      sdssr[count] = mr[i];
      ecatid[count] = sed_ids[i];
      dist8k[count] = dist8[i];
      nndistk[count] = nndist[i];
      nndist_percentk = nndist_percent[i];
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

  cout << "Number of galaxies with idx == true: " << count << endl;
  cout << "Number of galaxies with shapes: " << keep << endl;
  assert(count==keep);

  cout << "Creating TRUTH HDU" << endl;
  static string ttablename("TRUTH");
  Table* newTable;

  try{
    newTable = tFits->addTable(ttablename,keep,tcolName,tcolForm,tcolUnit);
  }
  catch(...){
    printf("Could not create table\n");
    exit(1);
  }

  cout << "Writing columns" << endl;
  try{
    cout << "writing id 1" << endl;
    newTable->column(tcolName[0]).write(id,1);
    cout << "writing id 2" << endl;
    newTable->column(tcolName[1]).write(id,1);
    cout << "writing sed_ids" << endl;
    newTable->column(tcolName[2]).write(ecatid,1);
    cout << "writing coeffs" << endl;
    newTable->column(tcolName[3]).write(cf,count,1);
    cout << "writing tmag" << endl;
    newTable->column(tcolName[4]).write(tmag,count,1);
    cout << "writing omag" << endl;
    newTable->column(tcolName[5]).write(omag,count,1);
    cout << "writing omag" << endl;
    newTable->column(tcolName[6]).write(flux,count,1);
    cout << "writing flux" << endl;
    newTable->column(tcolName[7]).write(ivar,count,1);
    cout << "writing ivar" << endl;
    newTable->column(tcolName[8]).write(omagerr,count,1);
    cout << "writing coeffs" << endl;
    newTable->column(tcolName[9]).write(am,count,1);
    cout << "writing ra" << endl;
    newTable->column(tcolName[10]).write(ra,1);
    cout << "writing dec" << endl;
    newTable->column(tcolName[11]).write(dec,1);
    cout << "writing z" << endl;
    newTable->column(tcolName[12]).write(z,1);
    cout << "writing haloid" << endl;
    newTable->column(tcolName[13]).write(haloid,1);
    cout << "writing 18" << endl;
    newTable->column(tcolName[18]).write(central,1);
    cout << "writing 21" << endl;
    newTable->column(tcolName[21]).write(e,count,1);
    cout << "writing 26" << endl;
    newTable->column(tcolName[26]).write(sdssr,1);
    cout << "writing 27" << endl;
    newTable->column(tcolName[27]).write(s,1);
    cout << "writing 28" << endl;
    newTable->column(tcolName[28]).write(px,1);
    cout << "writing 29" << endl;
    newTable->column(tcolName[29]).write(py,1);
    cout << "writing 30" << endl;
    newTable->column(tcolName[30]).write(pz,1);
    cout << "writing 31" << endl;
    newTable->column(tcolName[31]).write(vx,1);
    cout << "writing 32" << endl;
    newTable->column(tcolName[32]).write(vy,1);
    cout << "writing 33" << endl;
    newTable->column(tcolName[33]).write(vz,1);
    cout << "writing 34" << endl;
    newTable->column(tcolName[34]).write(e,count,1);
    cout << "writing 35" << endl;
    newTable->column(tcolName[35]).write(s,1);
    cout << "writing 38" << endl;
    newTable->column(tcolName[38]).write(dist8k,1);
    cout << "writing 39" << endl;
    newTable->column(tcolName[39]).write(nndistk,1);
    cout << "writing 40" << endl;
    newTable->column(tcolName[40]).write(nndist_percentk,1);
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
