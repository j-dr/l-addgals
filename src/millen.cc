#include <vector>
#include <cassert> 
#include <unistd.h>     //for sleep()
#include <fstream>
#include "hv2.h"
#include "hv.h"

#define COLORS
#define PRINTHALOS
//#define NOCFS


#ifdef HVL
Simulation sim(box, use_cells);
#endif
#ifdef LANL
Simulation sim(box);
#endif

Cosmology cosmo = sim.SimCosmology();
using namespace std;

extern "C" void covar(float *x1, float *y1, float *z1, 
		      float *vx1, float *vy1, float *vz1, int np1,
		      float rcube, float rmin, float rmax, int nbin);

void DeleteAndNullifyHighzgal(Galaxy*& pgalaxy){
  if((pgalaxy->Z()>ZREDMAX)||(pgalaxy->P()->Dec()>DECMAX)||(pgalaxy->P()->Ra()>RAMAX)
     ||(pgalaxy->Z()<ZREDMIN)||(pgalaxy->P()->Dec()<DECMIN)||(pgalaxy->P()->Ra()<RAMIN)
     ||(pgalaxy->P()->X()<=0)&&(pgalaxy->P()->Y()<=0)&&(pgalaxy->P()->Z()<=0)){
    delete pgalaxy;
    pgalaxy =0;
  }
}

int main(void){
  cosmo.Print();
  void AssignGalaxies(vector <Particle *> &particles, vector <Galaxy *> &galaxies);
  void UnAssignGalaxies(vector <Particle *> &particles, vector <Galaxy *> &galaxies);
  void SwapGalaxies(vector <Galaxy *> &galaxies, vector <Halo *> &halos);

  //make directory if it doesn't exist
  string mkstring = "mkdir "+out_path+"; cp ./box.h ./constants.h ./chainel.h ./chainstart.dat ./denspdf/hv_denspdf.dat "+out_path;
  PRNT("hv",mkstring);
  system(mkstring.c_str());

  //  ChainEl chainel();
  ChainEl chainel(-0.475, -0.2517, 2.856, 0.0, 2.2, -0.0536, 0.8785); 

  // Read in z as a function of distance so that redshifts 
  // can be determined without integrating
  cosmo.ReadZFile();
  ReadLFFile();

  double ngbar = LumNumberDensity(Magmin);

  vector <Galaxy *> galaxies;
  // Get a list of galaxies, whose magnitudes are determined by 
  // integrating the luminosity function
  // The number is determined by the constant "dim", 
  // which specifies the volume.
  // The galaxy constructor sets the properties:  cluster_gal, d8

  MSG("[hv] Getting galaxy luminosities and local densities");
  PRNTV(volume_fraction);
  PRNTV(sim.Boxsize());
  PRNTV(sim.CubeVolume());
  //  galaxies=GetGalaxies(sim.CubeVolume()*volume_fraction);

  
  galaxies=GetGalaxies(sim.CubeVolume()*volume_fraction, chainel);
  //  string tmpoutfile = out_path+"gdenstest.dat";
  //ofstream tmpout(tmpoutfile.c_str());
  //for(int gi=0;gi<galaxies.size();gi++){
  //if(gi!=galaxies[gi]->Gid())
  //  cout<<gi<<" "<<galaxies[gi]->Gid()<<endl;
  // }
  //tmpout<<galaxies[gi]->Dist8()<<" "<<galaxies[gi]->Mr()<<endl;
  
  //tmpoutfile = out_path+"pdenstest.dat";
  //ofstream tmpout2(tmpoutfle.c_str());
  //for(int pi=0;pi<galaxies.size();pi++)
  // tmpout2<<particles[pi]->Dist8()<<endl;
  
  cout<<"[hv] Assigning "<<galaxies.size()
      <<" galaxies to "<<particles.size()<<" particles ("
      <<galaxies.size()*1.0/particles.size()<<")"<<endl;

  AssignGalaxies(particles, galaxies);
  cout<<"swapping"<<endl;
  // This function also assigns gals to halos and calculates HOD.
  SwapGalaxies(galaxies, halos);

  // Remove galaxies that are outside of the z/dec/ra range.  
  PRNTVS("assigngals",galaxies.size(),"...removing high z/out of range galaxies");
  for_each(galaxies.begin(),galaxies.end(),DeleteAndNullifyHighzgal);
  galaxies.erase(remove(galaxies.begin(),galaxies.end(),static_cast<Galaxy*>(0)), galaxies.end());
  MSG(galaxies.size());

  int central_galaxies = 0;
  for(int gi=0;gi<galaxies.size();gi++){
    if(galaxies[gi]->Central())
      central_galaxies++;
  }
  cout<<"central galaxies in after removal:"<<central_galaxies<<endl;


  HaloOcc(halos);
  PRNTV("galaxies.size()");
  
    // Delete the unchosen particles to save memory.
    // But beware:  don't try to access these particles!!

    //NOTE: this SHOULD work.  but it seems to be introducing various bugs.
    //cout<<"Erasing unchosen particles"<<endl;
    //for_each(particles.begin(),particles.end(),DeleteAndNullifyUnchosenParticle);
    //particles.erase(remove(particles.begin(),particles.end(),static_cast<Particle*>(0)), particles.end());
   
  PRNT("hv",galaxies.size());
  ofstream outpfile(outpfn.c_str());
  ofstream outdfile(outdfn.c_str());
  ofstream outgfile(outgfn.c_str());
  ofstream outghfile(outghfn.c_str());
  ofstream outcfile(outcfn.c_str());
  
  //string tmpname = out_path+"z.dat";
  //  ofstream tmpfile(tmpname.c_str());

  cout<<"[hv] Writing cf file"<<endl;
  for(int gi=0;gi<galaxies.size();gi++){
    galaxies[gi]->WriteCFinfo(outcfile);
    //tmpfile<<galaxies[gi]->Z()<<" "<<galaxies[gi]->P()->R()<<endl;
  }

  //#ifndef COLORS
    string outstring = out_path+"dd.dat";
    ofstream outdd(outstring.c_str());
    for(int gi=0;gi<galaxies.size();gi++){
      if(randbool(0.1))
	outdd<<galaxies[gi]->Dist8()<<" "<<endl;
    }
    //#endif

  //Note:  new version assumes that nearest neighbor distances
  //are calculated for *all* galaxies in the data.
  cout<<"[hv] Getting SDSS galaxies."<<endl;
  system("date");
  //  vector <SEDTuple> SEDs;
  vector <GalSED> galseds = ReadSED();

  //#ifdef NEIGHBORS  
  cout<<"[hv] Getting nearest neighbors."<<endl;
  system("date");
  vector <Galaxy *> galaxycopy = galaxies;
  //  vector <float> nndist = GetNeighborDist2(galaxycopy);
  vector <float> nndist = GetNeighborDist(galaxycopy);
  // vector <float> nndist = GetNeighborDist3D(galaxycopy);
  string tmpoutfile = out_path+"dmdgmr.dat";
  ofstream outddmfile(tmpoutfile.c_str());
  cout<<galaxies.size()<<endl;
  for(int gi=0;gi<galaxies.size();gi++){
    Galaxy * gal = galaxies[gi];
    outddmfile<<gal->Dist8()<<" "<<nndist[gi]<<" "<<gal->Mr()<<endl;
  }

  //#ifdef COLORS

  cout<<"[hv] Assigning colors."<<endl;
  system("date");
  galaxycopy = galaxies;
  vector <int> sed_ids = GetSEDs(galaxycopy, nndist, galseds);
  cout<<"[hv] Printing "<<galaxies.size()<<" galaxies "<<endl;
  system("date");
  // Passivly evolve galaxy magnitudes
  MSG("Evolving galaxies now");
  //if(evolution == BLAN)
    for_each(galaxies.begin(),galaxies.end(),EvolveGal);
  system("date");
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

  MSG("[hv] Printing galaxies in volume");
  cout<<galaxies.size()<<endl;
  for(int gi=0;gi<galaxies.size();gi++){
    Galaxy * gal = galaxies[gi];
    Particle * p = gal->P();
    assert(p);  //this better be true since you removed the other ones.
    int hid = p->Hid();
    int sedid = sed_ids[gi];
    if((nndist[gi]>0)&&
       (p->X()>=0)&&(p->Y()>=0)&&(p->Z()>=0)&&(sedid>=0)){
      //      outgfile<<sed_ids[gi]<<" ";
      outgfile<<galseds[sed_ids[gi]].CatId()<<" ";
      //      SEDs[sed_ids[gi]].Write(outgfile);
      gal->Write(outgfile);
      p->Write(outpfile);
      if(hid<0){
	outghfile<<"0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "<<endl;
      }
      else{
	outghfile//<<in_vol<<" "
	  <<hid<<" "
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
      outdfile<<gal->Dist8()<<endl;
    }
    else{
      cout<<"[hv] didn't print galaxy "
	  <<gi<<" "
	  <<nndist[gi]<<" "
	  <<sedid<<" "
	  <<p->X()<<" "
	  <<p->Y()<<" "
	  <<p->Z()<<endl;
    }
  }
  MSG("Exiting normally");
}
