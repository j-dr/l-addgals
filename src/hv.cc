#include <vector>
#include <cassert> 
#include <unistd.h>     //for sleep()
#include <fstream>
#include <math.h> 
#include <cstring>
#include <stdlib.h>
#include "hv2.h"
#include "hv.h"
#include "kcorrect_utils.h"
#include "errormodel.h"
#include "shapes.h"
#include "fivetuple.h"
#include "ReadParameters.h"
#include "global_vars.h"
#include "outputs.h"
#include "galaxy_global.h"
#include "galaxy.h"
#ifdef HEALPIX
#include "healpix_utils.h"
#endif

//#define COLORS
//#define PRINTHALOS
//#define NOCFS

using namespace std;

extern "C" void covar(float *x1, float *y1, float *z1, 
		      float *vx1, float *vy1, float *vz1, int np1,
		      float rcube, float rmin, float rmax, int nbin);
extern double normal_random(float mean, float stddev);
extern double ranf(void);

void DeleteAndNullifyHighzgal(Galaxy*& pgalaxy){
#ifdef FULL_SKY
  if((pgalaxy->Z()>ZREDMAX)||(pgalaxy->P()->Dec()>DECMAX)||(pgalaxy->P()->Ra()>RAMAX)
     ||(pgalaxy->Z()<ZREDMIN)||(pgalaxy->P()->Dec()<DECMIN)||(pgalaxy->P()->Ra()<RAMIN)){
    delete pgalaxy;
    pgalaxy =0;
  }
#else
  if((pgalaxy->Z()>ZREDMAX)||(pgalaxy->P()->Dec()>DECMAX)||(pgalaxy->P()->Ra()>RAMAX)
     ||(pgalaxy->Z()<ZREDMIN)||(pgalaxy->P()->Dec()<DECMIN)||(pgalaxy->P()->Ra()<RAMIN)
     ||(pgalaxy->P()->X()<=0)&&(pgalaxy->P()->Y()<=0)&&(pgalaxy->P()->Z()<=0)){
    delete pgalaxy;
    pgalaxy =0;
  }
#endif
}

void DeleteAndNullifyHighzDimgal(Galaxy*& pgalaxy){
  if((pgalaxy->zGal()>ZREDMAX)||(pgalaxy->zGal()<ZREDMIN)){
    delete pgalaxy;
    pgalaxy =0;
  }
}

//Added by mbusha
void DeleteAndNullifyCentrals(Galaxy*& pgalaxy){
  if(pgalaxy->Central()){
    delete pgalaxy;
    pgalaxy =0;
  }
}
void add_bcgs(vector <Particle *> &particles, vector <Galaxy *> &galaxies, vector  <Halo *> &halos);
void AssignBCGs(vector <Particle *> &particles, vector <Galaxy *> &galaxies, vector <Halo *> &halos);
bool DLessGal(Galaxy * a, Galaxy * b);
bool MLessGal(Galaxy * a, Galaxy * b);
void MakeLL(int &LLBins, int * &NInBin, int * &BinStart, vector <Particle *> &particles);
int GalaxyZBin(float zRed);
void PrintMinMaxDens(vector <Galaxy *> &galaxies, vector <Particle *> &particles);
void print_galaxies(vector <Galaxy *> &galaxies, vector <Particle *> &particles, vector <Halo *> &halos, vector <GalSED> galseds, vector <int> sed_ids, vector <float> nndist, vector <float> nndist_percent, string outpfn, string outdfn, string outgfn, string outghfn, string outgzfn, string outrfn);
//void print_galaxies(vector <Galaxy *> &galaxies, vector <Particle *> &particles, vector <Halo *> &halos, vector <GalSED> galseds, vector <int> sed_ids, vector <float> nndist);

//added by mbusha to make sure the nn operation works correctly
class BrightGal{
public:
  BrightGal(float junk):_junk(junk){}
  bool operator()(Galaxy* g) const;
  
private:
  float _junk;
};

inline bool BrightGal::operator()(Galaxy* g) const {
  return (evolve_mag(g->Mr(),g->Z())<Magmin_dens);
}

#ifdef COLORS_FROM_RELATIVE_DENSITY
vector <float> GetNeighborPercents(vector <float> nndist, vector <Galaxy *> &galaxies);
#endif

#ifdef JUST_COLORS
void read_galaxies(vector <Particle *> &particles, vector <Galaxy *> &galaxies, vector <Halo *> &halos);
#endif
vector <Galaxy *> ReadSHAMGalaxies(string shamfile);
void get_mhost(vector <Halo *> &halos, vector <Galaxy *> &galaxies);

void set_d8(vector <Galaxy *> &galaxies, int idenspdf);
void set_d8(vector <Galaxy *> &galaxies, FiveTuple ftup);

Simulation DefineSimulation(void);

#ifdef HEALPIX
//Calculates the area in a halpix pixel and our ra/dec limits by
//just laying down a bunch of random points on the sky and measuring
//the fraction that are in the specified limits and pixel.  
float healpix_sqdeg(void)
{
  long npoints = 100000000;
  long n_in_area = 0;
  long ThisPixel = 0;
  int i;
  float theta, phi, theta_dec, fraction;

  for(i=0;i<npoints;i++)
    {
      phi = 360.0*drand48();
      if (phi < RAMIN || phi > RAMAX)
	continue;
      theta = acos(drand48()*2.0 - 1.0);
      theta_dec = 90.0 - (theta*180.0/PI);
      if(theta_dec < DECMIN || theta_dec > DECMAX)
	continue;
      phi *= (PI/180.0);
      ang2pix_ring(nSide, theta, phi, &ThisPixel);
      if(ThisPixel == PixelNum)
	n_in_area++;
    }
  fraction = ((float) n_in_area)/((float) npoints);
  //fraction /= 10.0; //for testing purposes we make things smaller/faster
  return fraction*41253.0;
}
#endif

//---------------------------------------------------------------------------\\
//          Main Code starts here                                            \\
//---------------------------------------------------------------------------\\

int main(void){
  double t1, t2, TimeAssign;

  cout<<endl<<"Thank you for using ADDGALS!  The code is now running..."<<endl<<endl;

#ifdef HEALPIX
  cout<<"Note:  Using healpix decomposition within ra/dec limits.  Besure nSide and PixelNum are defined in the NumericalParameters file."<<endl<<endl;
  //float kara = 301.685;
  //float kadec = 4.9077;
  //long kapixel;
#endif
#ifdef RED_FRACTION
  cout<<"Assuming red fraction evolution model..."<<endl;
#endif

  //First we read in our parameters
  readParameters();
  //ang2pix_ring(nSide, (90.-kadec)*3.1415926/180.0, kara*3.1415926/180.0, &kapixel);
  //cout<<kapixel<<endl;
  cout<<"Status of read_hod: "<<read_hod<<endl;
#ifdef HEALPIX
  cout<<"Specified PixelNum: "<<PixelNum<<endl;
#endif

  cout<<"Using evolution model: "<<evolution;
  if (evolution == 0) cout<<" (NOEV)"<<endl;
  if (evolution == 1) cout<<" (BLAN)"<<endl;
  if (evolution == 2) cout<<" (FABER)"<<endl;
  if (evolution == 3) cout<<" (TIME)"<<endl;

  cout << "ring2nest test " << ring2nest(PixelNum, 2) << endl;

  //correct out output files
  outpfn = out_path+outpfn;
  outdfn = out_path+outdfn;
  outgfn = out_path+"/"+flabel+"."+PSTR+"."+ZSTR+".fits";;
  outghfn= out_path+outghfn;
  outhfn = out_path+outhfn;
  outcfn = out_path+outcfn;
  outgzfn = out_path+outgzfn;
  outrfn = out_path+outrfn;
  outafn = out_path+outafn;

  cout << "outgfn: " << outgfn << endl;

  sim = DefineSimulation();
  cosmo = sim.SimCosmology();

  //Print some basic information on the cosmology and type of simulation being used
  cosmo.Print();
  cout<<"Simulation definition:  "<<sim.Type()<<endl;
  PRNTV(sim.Boxsize());

  // Tabulate the redshift-distance relation
  //cosmo.GetZofR(sim.OmegaM(), sim.OmegaL());
  //do everything assuming a fiducial cosmological model
  cosmo.GetZofR(0.3, 0.7);

  //Calculate the volume of the region.  Convert to radians and switch dec to standard spherical units
  float r_zmin = cosmo.RofZ(ZREDMIN);
  float r_zmax = cosmo.RofZ(ZREDMAX);
#ifdef HEALPIX
  //need to calculate the area in a pixel cell
  float area = healpix_sqdeg(); 
  cout<<"  Sq deg in this pixel: "<<area<<endl;
  float volume = 4./3*PI*(pow(r_zmax,3)-pow(r_zmin,3))*area/41253.;
#else
  float decmin_rad = M_PI*(90 - DECMAX)/180;
  float decmax_rad = M_PI*(90 - DECMIN)/180;
  float ramin_rad = M_PI*RAMIN/180;
  float ramax_rad = M_PI*RAMAX/180;
  float dra = ramax_rad - ramin_rad;
  float volume = (ramax_rad-ramin_rad)*(cos(decmin_rad)-cos(decmax_rad))*(pow(r_zmax,3)-pow(r_zmin,3))/3.;
#endif
#ifdef SNAPSHOT
  volume = sim.Boxsize()*sim.Boxsize()*sim.Boxsize();
#endif
  RMIN_REAL = r_zmin;
  RMAX_REAL = r_zmax;
  cout<<"Specified Simulation Volume: "<<volume<<endl;
  cout<<"Radial limits: "<<RMIN_REAL<<" - "<<RMAX_REAL<<endl;
  cout<<"Redshift limits: "<<ZREDMIN<<" "<<ZREDMAX<<endl;


  //Define some functions that are called later in the code
#ifdef SHAM_TEST
  void Assignment(vector <Particle *> &particles, vector <Galaxy *> &galaxies, vector <Halo *> &halos);
#else
  void Assignment(vector <Particle *> &particles, vector <Galaxy *> &galaxies);
#endif
  void UnAssignGalaxies(vector <Particle *> &particles, vector <Galaxy *> &galaxies);
  void SwapGalaxies(vector <Galaxy *> &galaxies, vector <Halo *> &halos);

  //declare our big data arrays -- they're empty for now.  
  vector <Particle*> particles;
  vector <Galaxy *> galaxies;
  vector <Halo*> halos;

  //make output directory if it doesn't exist
  string mkstring = "mkdir "+out_path+"; cp ./NumericalParameters ./StringParameters "+out_path;
  PRNT("hv",mkstring);
  system(mkstring.c_str());

#ifdef LF_FROM_DATA
  read_lf_data();
#endif

  //ReadLFFile(); //pretty sure this does nothing

  cout<<endl;

#ifndef JUST_COLORS  //We can skip to color assignment if we read in full galaxy list

  //Calculate some basic properties about our galaxy distribution
  double ngbar = LumNumberDensity(Magmin);
  PRNT("hv",LumNumberDensity(Magmin));
  //mbusha:  I don't think LumNumberDensityInterp is ever called, and it keeps failing for some reason!
  //PRNT("hv",LumNumberDensityInterp(Magmin));
  PRNT("hv",Magmin);

#ifndef TESTREADPART
  //Read the halos and particles from the input files
  MSG("[hv] Reading halos");
  halos = ReadHalos();
#endif

  MSG("[hv] Reading particles");
  particles = ReadParticles();
  PRNTV(particles.size());

  cout << "PARTICLE_FRACTION: " << PARTICLE_FRACTION << " Particle Mass: " << sim.ParticleMass() << " Volume: " << volume << endl;
  float rho = (particles.size()/PARTICLE_FRACTION)*sim.ParticleMass()/volume;
  float rhoback = 3.*100.0*100.0/(8*M_PI*4.301e-9)*sim.OmegaM();
  float phi_rescale = rho/rhoback;
  cout<<"rho and rhoback for cosmic variance calculaiton: "<<rho<<" "<<rhoback<<endl;
  cout<<"Rescaling factor for LF normalization: "<<phi_rescale<<endl;


#ifdef SNAPSHOT
  phi_rescale = 1.0;
  float posmin = 1000.;
  float posmax = -1000.;
  for(int ip=0;ip<particles.size();ip++)
  {
    if (particles[ip]->X() < posmin) posmin = particles[ip]->X();
    if (particles[ip]->Y() < posmin) posmin = particles[ip]->Y();
    if (particles[ip]->Z() < posmin) posmin = particles[ip]->Z();
    if (particles[ip]->X() > posmax) posmax = particles[ip]->X();
    if (particles[ip]->Y() > posmax) posmax = particles[ip]->Y();
    if (particles[ip]->Z() > posmax) posmax = particles[ip]->Z();
  }
  float boxsize = posmax - posmin;
  volume = boxsize*boxsize*boxsize;
  cout<<"Boxsize and volume as calculated from particle distribution: "<<boxsize<<", "<<volume<<endl;
#endif

  cout<<"Particles read: check mem"<<endl;
  system("ps ux | grep hv > mem.tmp");
  cout<<"[hv]"<<particles.size()<<" "<<particles.max_size()<<" "<<particles.capacity()<<endl;
  cout<<"[hv] BOXSIZE:"<<sim.Boxsize()<<endl;

  MSG("Particles Read In: check mem");
  cout<<"[hv] size:"<<particles.size()<<" max:"<<particles.max_size()<<" capacity:"<<particles.capacity()<<endl;
  system("ps ux | grep hv >> mem.tmp");

  //float rho = particles.size()*sim.ParticleMass()/volume;
  //float rhoback = 3.*100*100/(8*!PI*4.301e-9)*sim.OmegaM();
  //float phi_rescale = rho/rhoback;
 
  double fraction_of_particles = ngbar*volume/particles.size();
  //PRNT("hv",fraction_of_particles);
  cout<<"Fraction of particles that will be assigned a galaxy: "<<fraction_of_particles<<endl;

  //Now we're done loading the particles



  // Sort particles on redshift.  Only do this if we're running a LC.   
  // mbusha:  Note that I don't think this is evern necessary anymore because a sort is done in assignment....
#ifndef SNAPSHOT
  //MSG("[hv] Sorting particles.");
  //sort(particles.begin(),particles.end(),zLess);
  //MSG("[hv] Particles sorted: check mem");
  //cout<<"[hv]"<<particles.size()<<" "<<particles.max_size()<<" "<<particles.capacity()<<endl;
  //cout<<endl;
#endif

  //Output some diagnostics on the rnn and redshift distribution of particles.
  //	This is just sanity checking.  
  minrnn = 10000.;
  maxrnn = 0.;
  float avgrnn = 0.;
  float maxz = 0.;
  float minz = 1000.;
  for(int i=0;i<particles.size();i++){
    if(particles[i]->Dist8() > maxrnn)
      maxrnn = particles[i]->Dist8();
    if(particles[i]->Dist8() < minrnn)
      minrnn = particles[i]->Dist8();
    if(particles[i]->Dist8() <= 0. || particles[i]->Dist8() > 25.0){
      cout<<"Error with particle "<<i<<"!  rnn = "<<particles[i]->Dist8()<<", ";
      particles[i]->PosPrint();
    }
    if(particles[i]->Zred() < minz) minz = particles[i]->Zred();
    if(particles[i]->Zred() > maxz) maxz = particles[i]->Zred();
    avgrnn += particles[i]->Dist8();
  }
  avgrnn /= particles.size();
  cout<<" minimum rnn = "<<minrnn<<", max rnn = "<<maxrnn<<", avg rnn = "<<avgrnn<<endl;
  cout<<" minimum z = "<<minz<<", maximum z = "<<maxz<<endl;


  // Get a list of galaxies, whose magnitudes are determined by 
  // integrating the luminosity function
  // The number is determined by the constant "dim", 
  // which specifies the volume.
  // The galaxy constructor sets the properties:  cluster_gal, d8
  MSG("[hv] Getting galaxy luminosities and local densities");
#ifdef SHAM_TEST
  galaxies=ReadSHAMGalaxies(sham_file);
  cout<<"Reading in the HOD information file."<<endl;
  if (read_hod) get_mhost(halos, galaxies);
  //cout<<"Some mhost values: "<<galaxies[0]->Mhost()<<" "<<galaxies[2]->Mhost()<<endl;
  /*
  //If there are more galaxies than particles we cut the file down to a filling fraction of 1.0
  if (galaxies.size() > particles.size())
    {
      cout<<"More galaxies than particles: "<<galaxies.size()<<".  Deleting the dimmest ones."<<endl;
      sort(galaxies.begin(),galaxies.end(),MLessGal);
      cout<<"Galaxies sorted: "<<galaxies[0]->Mr()<<", "<<galaxies[galaxies.size()-1]->Mr()<<endl;
      galaxies.erase(galaxies.begin()+particles.size()*0.5,galaxies.end());
      cout<<"Number of remain galaxies: "<<galaxies.size()<<endl;
      cout<<" Magnitude range: "<<galaxies[0]->Mr()<<", "<<galaxies[galaxies.size()-1]->Mr()<<endl;
    }
  */
#elif OLDPDF
  galaxies=GetGalaxies(volume);
#else
  galaxies=GetGalaxies(volume, phi_rescale);
#endif

#ifdef BCGS
  //Add our BCGs to the galaxy distribution
  if (halos.size() > 0 && read_hod == 0) AssignBCGs(particles, galaxies, halos);
#endif

  //Sort galaxies by density for rapid assignment
  sort(galaxies.begin(),galaxies.end(),DLessGal);
  
  //do the galaxy assignmant
  cout<<"[hv] Assigning "<<galaxies.size()
      <<" galaxies to "<<particles.size()<<" particles (Filling Fraction: "
      <<galaxies.size()*1.0/particles.size()<<")"<<endl;
  PrintMinMaxDens(galaxies, particles);

  cout<<"Assigning galaxies"<<endl;
  t1 = clock();
#ifdef SHAM_TEST
  cout<<"read_hod = "<<read_hod<<endl;
  Assignment(particles, galaxies, halos);
#else
  Assignment(particles, galaxies);
#endif

#ifdef DEBUG
  //save properties of the galaxy and particle densities for later diagnostics
  ofstream outdfile(outdfn.c_str());
  for(int ig=0;ig<galaxies.size();ig++)
    outdfile<<galaxies[ig]->Dist8()<<" "<<galaxies[ig]->P()->Dist8()<<endl;

  string outstring = out_path+"dd.dat";
  ofstream outdd(outstring.c_str());
  for(int gi=0;gi<galaxies.size();gi++){
    if(randbool(0.1))
	outdd<<galaxies[gi]->Dist8()<<" "<<endl;
  }

#endif

  cout<<endl<<"Done with the assigning."<<endl;
  t2 = clock();
  TimeAssign = (t2-t1)/CLOCKS_PER_SEC;
  cout<<"Assigned galaxies in "<<TimeAssign<<" seconds."<<endl;
  cout<<endl;

  //This function doesn't actually do anything anymore -- but I might
  //want to go back and revisit it while adding the CLF stuff
#ifdef DEBUG
  cout<<"swapping"<<endl;
  // This function also assigns gals to halos and calculates HOD.
  SwapGalaxies(galaxies, halos);
  cout<<"Did the swap"<<endl;
#endif

  // Remove galaxies that are outside of the z/dec/ra range.  
  PRNTVS("hv",galaxies.size(),"...removing high z/out of range galaxies");
  for_each(galaxies.begin(),galaxies.end(),DeleteAndNullifyHighzgal);
  galaxies.erase(remove(galaxies.begin(),galaxies.end(),static_cast<Galaxy*>(0)), galaxies.end());
  MSG(galaxies.size());

  int central_galaxies = 0;
  for(int gi=0;gi<galaxies.size();gi++){
    if(galaxies[gi]->Central())
      central_galaxies++;
  }
  cout<<"central galaxies in after removal:"<<central_galaxies<<endl;

  //I dont' think this subroutine does anything important, but should confirm
  //HaloOcc(halos);
   
  PRNT("hv",galaxies.size());

#else //skip the if not JUST_COLORS statement
  read_galaxies(particles, galaxies, halos);
#endif

  //Read in our SED training set
  cout<<"[hv] Getting SDSS galaxies."<<endl;
  system("date");
  t1 = clock();
  vector <GalSED> galseds = ReadSED();
  vector <int> sed_ids;

#ifndef COLORS
  //If we're not doing colors, just make some dummy arrays
  vector <float> nndist(galaxies.size());
  vector <float> nndist_percent(galaxies.size());
  for(int i=0;i<galaxies.size();i++){
    nndist[i] = 0.;
    nndist_percent[i] = 0.;
  }
#else
  //get the colors
  cout<<"[hv] Getting nearest neighbors."<<endl;
  system("date");

  //this version just uses galaxies brighter than omagmin_dens
  cout<<" Making array of bright galaxies."<<endl;
  float junk = 0.;
  BrightGal bg(junk);
  cout<<" Doing any magnitude cuts"<<endl;
  vector <Galaxy *> galaxycopy;
  cout<<" doing copy_if"<<endl;
  for(int i=0;i<galaxies.size();i++)
    if (galaxies[i]->Mr() <= Magmin_dens)
      galaxycopy.push_back(galaxies[i]);


  //Calculate the projected distance
  cout<<"Getting Neighbors..."<<endl;
  vector <float> nndist = GetNeighborDist(galaxycopy, galaxies);

  cout<<"Got neighbor dists: check mem"<<endl;
  system("ps ux | grep hv > mem.tmp");

#ifdef  COLORS_FROM_RELATIVE_DENSITY
  cout<<"Getting Neighbor Percents..."<<endl;
  vector <float> nndist_percent = GetNeighborPercents(nndist, galaxies);
  cout<<"Measured the nndist_percent's: "<<nndist_percent[0]<<" "<<nndist_percent[1]<<endl;
  cout << "check mem" << endl;
  system("ps ux | grep hv > mem.tmp");
  //  ofstream outrfile(outrfn.c_str());
#ifdef DEBUG
  ofstream outrfile(outrfn.c_str());
  for(int i=0;i<galaxies.size();i++)
    outrfile<<nndist[i]<<" "<<nndist_percent[i]<<endl;
#endif
#else
  vector <float> nndist_percent(galaxies.size());
#endif

  cout<<"Have Neighbors"<<endl;
  float min_nndist = 1e10;
  float max_nndist = 0.0;
  for (int i=0;i<galaxies.size();i++){
    if (nndist[i] > max_nndist) 
      max_nndist = nndist[i];
    if (nndist[i] < min_nndist) 
      min_nndist = nndist[i];
  }
  cout<<"min/max nndist = "<<min_nndist<<"/"<<max_nndist<<endl;

  int Nnndist0 = 0;
  for(int inn=0;inn<nndist.size();inn++)
    if(nndist[inn] == 0.)
      Nnndist0++;
  cout<<"Number of galaxies with a nndist of 0:  "<<Nnndist0<<endl;
  string tmpoutfile = out_path+"dmdgmr.dat";
  ofstream outddmfile(tmpoutfile.c_str());
  cout<<galaxies.size()<<endl;
  for(int gi=0;gi<galaxies.size();gi++){
    Galaxy * gal = galaxies[gi];
    outddmfile<<gal->Dist8()<<" "<<nndist[gi]<<" "<<gal->Mr()<<endl;
  }

  galaxycopy = galaxies;
  cout<<"[hv] Assigning colors."<<endl;
  system("date");
  cout << "check mem" << endl;
  system("ps ux | grep hv > mem.tmp");

#ifndef COLORS_FROM_RELATIVE_DENSITY
  sed_ids = GetSEDs(galaxycopy, nndist, galseds, halos);
#else
  sed_ids = GetSEDs(galaxycopy, nndist_percent, galseds, halos);
#endif

  cout<<"[hv] Printing "<<galaxies.size()<<" galaxies "<<endl;
  system("date");

#endif // doing the colors?

  //Passivly evolve galaxy magnitudes
  MSG("Evolving galaxies now");
  for_each(galaxies.begin(),galaxies.end(),EvolveGal);
  system("date");

  t2 = clock();
  double TimeColor = (t2-t1)/CLOCKS_PER_SEC;
  cout<<"Assigned galaxy SEDs in "<<TimeColor<<" seconds."<<endl;

#ifdef PRINTHALOS
  MSG("[hv] Printing halos in volume");
  ofstream outhfile(outhfn.c_str());
  for(int hi=0;hi<halos.size();hi++){
    Halo * h = halos[hi];
    if(h->InVol())
      outhfile<<h->Id()<<" "
	      <<h->M()<<" "<<halos[hi]->R200()<<" "
	      <<h->Ra()<<" "<<h->Dec()<<" "<<h->Zred()<<" "<<h->Ngal()<<endl;
  }
#endif
#ifdef BCC
  //Eventually need to add u
  int nbands = 5;
  int nelem = 3;
  int ntemp = 5;
  float band_shift = 0.1;
  char filterfile[2000];
  vector<float> mr(galaxies.size());
  vector<float> z(galaxies.size());
  vector<int> id(galaxies.size());
  vector<float> coeff(galaxies.size()*ntemp);
  vector<float> tmag(galaxies.size()*nbands);
  vector<float> amag(galaxies.size()*nbands);
  read_out_galaxy_info(galaxies, mr, z, galseds, sed_ids, id);
  match_coeff(id, &coeff[0]);

  strcpy(filterfile, "/nfs/slac/g/ki/ki23/des/jderose/l-addgals/src/des_filters.txt");

  float zmin = 0.0;
  float zmax = 2.5;

  assign_colors(mr, coeff, z, zmin, zmax,
		band_shift, nbands, filterfile, tmag, amag);

  vector<float> omag(galaxies.size()*nbands);
  vector<float> omagerr(galaxies.size()*nbands);
  vector<float> flux(galaxies.size()*nbands);
  vector<float> ivar(galaxies.size()*nbands);
  vector<bool> idx(galaxies.size(), true);

  observe_des_y5(tmag, flux, ivar, omag, omagerr, idx);
  assert(omag.size()%nbands == 0);
  vector<double> e(omag.size()/nbands*2);
  vector<double> s(omag.size()/nbands);
  generate_shapes(omag, e, s, nelem, nbands);
  
  write_bcc_catalogs(galaxies, particles, amag, tmag, mr, 
		     omag, omagerr, flux, ivar, e, s,
		     idx, halos, id, coeff, outgfn, 
		     outghfn);
#else
  MSG("[hv] Printing galaxies in volume");
  cout<<galaxies.size()<<endl;
  print_galaxies(galaxies, particles, halos, galseds, sed_ids, nndist, nndist_percent, outpfn, outdfn, outgfn, outghfn, outgzfn, outrfn);

#endif

  MSG("Exiting normally");
}
