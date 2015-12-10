#include <vector>
#include <cassert> 
#include <unistd.h>     //for sleep()
#include <fstream>
#include "galaxy.h"

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

void  write_bcc_catalogs(vector<Galaxy *> &galaxies, vector<Particle *> &particles, 
			 vector<float> amag, vector<float> tmag, vector<float> mr, 
			 vector<float> omag, vector<float> omagerr, vector<float> flux, 
			 vector<float> fluxerr, vector<float> e, vector<float> s,
			 vector<bool> idx, vector <Halo *> &halos,
			 vector<int> &sed_ids, vector<float> &coeffs,
			 string outgfn, string outghfn)
{
  using namespace CCfits;
  std::unique_ptr<FITS> tFits;
  int size = galaxies.size();
  int keep = s.size()
  
  //Write truth file
  try{
    tFits.reset(new FITS(outgfn,Write));
  }
  catch (CCfits::FITS::CantOpen){
    cerr << "Can't open " << outfile << endl;
    opened = false;       
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

  tcolUnit[0] = "";
  tcolUnit[1] = "";
  tcolUnit[2] = "";
  tcolUnit[3] = "";
  tcolUnit[4] = "mag";
  tcolUnit[5] = "mag";
  tcolUnit[6] = "nmgy";
  tcolUnit[7] = "nmgy^2";
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
  
  tcolForm[0] = "K"
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

  vector<float> ra(keep), dec(keep), px(keep), py(keep), pz(keep),
    vx(keep), vy(keep), vz(keep), sdssr(keep), z(keep), id(keep);
  
  //Go through the galaxies and get the info we want
  int count=0;
  for (int i=0; i<size; i++)
    {
      if (!idx[i]) continue;

      ra[count] = galaxies[i]->Ra();
      dec[count] = galaxies[i]->Dec();
      id[count] = galaxies[i]->particle->Gid();
      px[count] = galaxies[i]->particle->X();
      py[count] = galaxies[i]->particle->Y();      
      pz[count] = galaxies[i]->particle->Z();
      vx[count] = galaxies[i]->particle->VX();
      vy[count] = galaxies[i]->particle->VY();
      vz[count] = galaxies[i]->particle->VZ();
      z[count] = galaxies[i]->Z();
      central[count] = galaxies[i]->Central()
      sdssr[count] = mr[i];

    }
  static string ttablename("TRUTH");
  Table* newTable;

  try{
    newTable = pFits->addTable(ttablename,keep,tcolName,tcolForm,tcolUnit);
  }
  catch(...){
    printf("Could not create table\n");
    exit(1);
  }

  try{
    newTable->column(colname[0]).write(id,1);
    newTable->column(colname[1]).write(id,1); //should this be something else?
    newTable->column(colname[2]).write(sed_id,1);
    newTable->column(colname[3]).write(coeff,1);
    newTable->column(colname[4]).write(tmag,1);
    newTable->column(colname[5]).write(omag,1);
    newTable->column(colname[4]).write(flux,1);
    newTable->column(colname[5]).write(ivar,1);
    newTable->column(colname[4]).write(omagerr,1);
    newTable->column(colname[5]).write(amag,1);
    newTable->column(colname[4]).write(ra,1);
    newTable->column(colname[5]).write(dec,1);
    newTable->column(colname[4]).write(z,1);
    newTable->column(colname[5]).write(haloid,1);
    newTable->column(colname[4]).write(rhalo,1);
    newTable->column(colname[5]).write(m200,1);
    newTable->column(colname[4]).write(ngals,1);
    newTable->column(colname[5]).write(r200,1);
    newTable->column(colname[4]).write(central,1);
    newTable->column(colname[5]).write(ra,1);
    newTable->column(colname[4]).write(dec,1);
    newTable->column(colname[5]).write(e,1);
    newTable->column(colname[4]).write(gamma1,1);
    newTable->column(colname[5]).write(gamma2,1);
    newTable->column(colname[4]).write(kappa,1);
    newTable->column(colname[5]).write(mu,1);
    newTable->column(colname[4]).write(mr,1);
    newTable->column(colname[5]).write(s,1);
    newTable->column(colname[4]).write(dec,1);
    newTable->column(colname[5]).write(e,1);
    newTable->column(colname[4]).write(px,1);
    newTable->column(colname[5]).write(py,1);
    newTable->column(colname[4]).write(pz,1);
    newTable->column(colname[5]).write(vx,1);
    newTable->column(colname[4]).write(vy,1);
    newTable->column(colname[5]).write(vz,1);
    newTable->column(colname[4]).write(e,1);
    newTable->column(colname[5]).write(s,1);
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

  //Write the halo information

}
