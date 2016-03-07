#include <vector>
#include <cassert> 
#include <memory>
#include <unistd.h>     //for sleep()
#include <fstream>
#include <CCfits/CCfits>
#include "outputs.h"

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
{
  return;
}

void write_bcc_catalogs(vector<float> &am, vector<float> &tmag, vector<float> &sdssr, 
			 vector<float> &omag, vector<float> &omagerr, vector<float> &flux, 
			 vector<float> &ivar, vector<double> &e, vector<double> &s,
			 vector<float> &ra, vector<float> &dec, vector<float> &px,
			 vector<float> &py, vector<float> &pz, vector<float> &vx,
			 vector<float> &vy, vector<float> &vz, vector<float> &z,
			 vector<int> &gid, vector<int> &central, vector<int> &haloid,
			 vector<int> &ecatid, vector<float> &cf,
			 string outgfn, string outghfn)
{
  using namespace CCfits;
  std::auto_ptr<FITS> tFits;
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

  vector<float> empty(keep, -1);

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
    newTable->column(tcolName[0]).write(gid,1);
    cout << "writing id 2" << endl;
    newTable->column(tcolName[1]).write(gid,1); //should this be something else?
    cout << "writing sed_ids" << endl;
    newTable->column(tcolName[2]).write(ecatid,1);
    cout << "writing coeffs" << endl;
    newTable->column(tcolName[3]).write(cf,keep,1);
    cout << "writing tmag" << endl;
    newTable->column(tcolName[4]).write(tmag,keep,1);
    cout << "writing omag" << endl;
    newTable->column(tcolName[5]).write(omag,keep,1);
    cout << "writing omag" << endl;
    newTable->column(tcolName[6]).write(flux,keep,1);
    cout << "writing flux" << endl;    
    newTable->column(tcolName[7]).write(ivar,keep,1);
    cout << "writing ivar" << endl;
    newTable->column(tcolName[8]).write(omagerr,keep,1);
    cout << "writing coeffs" << endl;
    newTable->column(tcolName[9]).write(am,keep,1);
    cout << "writing ra" << endl;
    newTable->column(tcolName[10]).write(ra,1);
    cout << "writing dec" << endl;
    newTable->column(tcolName[11]).write(dec,1);
    cout << "writing z" << endl;
    newTable->column(tcolName[12]).write(z,1);
    cout << "writing haloid" << endl;
    newTable->column(tcolName[13]).write(haloid,1);
    cout << "writing 14" << endl;
    newTable->column(tcolName[14]).write(empty,1);
    cout << "writing 15" << endl;
    newTable->column(tcolName[15]).write(empty,1);
    cout << "writing 16" << endl;
    newTable->column(tcolName[16]).write(empty,1);
    cout << "writing 17" << endl;
    newTable->column(tcolName[17]).write(empty,1);
    cout << "writing 18" << endl;
    newTable->column(tcolName[18]).write(central,1);
    cout << "writing 19" << endl;
    newTable->column(tcolName[19]).write(ra,1);
    cout << "writing 20" << endl;
    newTable->column(tcolName[20]).write(dec,1);
    cout << "writing 21" << endl;
    newTable->column(tcolName[21]).write(e,keep,1);
    cout << "writing 22" << endl;
    newTable->column(tcolName[22]).write(empty,1);
    cout << "writing 23" << endl;
    newTable->column(tcolName[23]).write(empty,1);
    cout << "writing 24" << endl;
    newTable->column(tcolName[24]).write(empty,1);
    cout << "writing 25" << endl;
    newTable->column(tcolName[25]).write(empty,1);
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
    newTable->column(tcolName[34]).write(e,keep,1);
    cout << "writing 35" << endl;
    newTable->column(tcolName[35]).write(s,1);
    cout << "writing 36" << endl;
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

int extract_galaxy_information(vector<Galaxy *> &galaxies, vector<Particle *> &particles, 
			       vector<float> &amag, vector<float> &mr, vector<bool> &idx, 
			       vector <Halo *> &halos, vector<float> &coeffs, 
			       vector<float> &ra, vector<float> &dec, vector<float> &px, 
			       vector<float> &py, vector<float> &pz, vector<float> &vx, 
			       vector<float> &vy, vector<float> &vz, vector<float> &sdssr, 
			       vector<float> &z, vector<int> &gid, vector<int> &central, 
			       vector<int> &haloid, vector<float> &cf, 
			       vector<float> &am, vector<int> &ecatid, int nbands,
			       int ntemp)
{
  cout << "Extracting relevant galaxy information" << endl;
  int count=0;
  int size=galaxies.size();
  Particle * p;
  for (int i=0; i<size; i++)
    {
      if (idx[i]) {
	p = galaxies[0]->P();
	int hid = p->Hid();
	ra.push_back( galaxies[0]->Ra() );
	dec.push_back( galaxies[0]->Dec() );
	gid.push_back( p->Gid() );
	px.push_back( p->X() );
	py.push_back( p->Y() );
	pz.push_back( p->Z() );
	vx.push_back( p->Vx() );
	vy.push_back( p->Vy() );
	vz.push_back( p->Vz() );
	central.push_back( galaxies[0]->Central() ) ;
	sdssr.push_back( mr[0] );
	ecatid[count] = ecatid[i];
	z[count] = z[i];
	for (int c=0; c<nbands; c++)
	  {
	    am.push_back( amag[c] );
	  }
	for (int c=0; c<ntemp; c++)
	  {
	    cf.push_back( coeffs[c] );
	  }
	if ((central[count]==1) & (hid>=0))
	  {
	    haloid.push_back( halos[hid]->Id() );
	  }
	else
	  {
	    haloid.push_back(-1);
	  }
	count++;
      }
      galaxies.erase(galaxies.begin());
      mr.erase(mr.begin());
      coeffs.erase(coeffs.begin(), coeffs.begin()+ntemp);
      amag.erase(amag.begin(), amag.begin()+nbands);
    }
  z.resize(count);
  ecatid.resize(count);
  cout << "Got information" << endl;
  return count;
}
 

