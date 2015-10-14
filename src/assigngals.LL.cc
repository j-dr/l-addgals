#include <vector>
#include <cassert> 
#include <iostream>
#include "galaxy.h"
#include "choose.h"
#include "point.h"
#include "particle.h"
#include "myrand.h"
#include "fivetuple.h"
#include "constants.h"
#include "halo.h"
#include "stl_util.h"  //for copy if


// this is when the galaxy didn't find a good particle.
// if this happens often you should worry
void DeleteAndNullifyBadgal(Galaxy*& pgalaxy)
{
  Particle * myp = pgalaxy->P();
  //  if((!myp)&&(!myp->Keep())){
  if((!myp)&&(pgalaxy)){
    cerr<<"removing unchosen galaxy: "
	<<pgalaxy->Mr()<<" "<<pgalaxy->Dist8()<<endl;
    delete pgalaxy;
    pgalaxy=0;
  }
}

// compare galaxy pointers by local mass density
bool DLessGal(Galaxy * a, Galaxy * b)
{
  return a->Dist8() < b->Dist8(); 
}

// compare galaxy pointers by local mass density
bool DLessP(Particle * a, Particle * b)
{
  return a->Dist8() < b->Dist8(); 
}

// unassigns galaxies from particlesand deletes galaxies
void UnAssignGalaxies(vector <Particle *> &particles, vector <Galaxy *> &galaxies){
  for(int i=0;i<particles.size();i++){
    particles[i]->UnMakeGal();
  }

  for(int i=0;i<galaxies.size();i++){
    delete galaxies[i];
    galaxies[i] = 0;
  }
  galaxies.erase(remove(galaxies.begin(),galaxies.end(),static_cast<Galaxy*>(0)), galaxies.end());
}

//void AssignGalaxies(vector <Particle *> &particles, vector <Galaxy *> &galaxies, vector <Halo *> &halos)

int GalaxyZBin(float zRed)
{
  float ThisZ = 0.0;
  int zBin = 0;
  while(ThisZ < zRed)
    {
      ThisZ += zTol*(1+ThisZ);
      zBin++;
    }
  zBin--;
  return zBin;
}

#ifdef 0
void MakeLL(int LLBins, int * &LLHoc, int * &LL, vector <Particle *> &particles)
{
  //Determine the numbeer of bins
  LLBins = GalaxyZBin(ZREDMAX)+1;
  cout<<"    The list will have "<<LLBins<<" bins."<<endl;

  //clear the linked list and HOC
  LLHoc = new int[LLBins];
  for(int i=0;i<LLBins;i++)
    LLHoc[i] = -1;
  for(int i=0;i<particles.size();i++)
    LL[i] = -1;

  float Percent = 0.0;
  //Add galaxies in z bins to linked list in order of increasing density
  for(int ip=0;ip<particles.size();ip++){
    if(float(ip)/float(particles.size()) > Percent){
      cout<<"    "<<100.*Percent<<"% done."<<endl;
      Percent += 0.01;
    }
    int ThisBin = GalaxyZBin(particles[ip]->Zred());
    //loop through the chain until we find the right Density to insert the particle
    int ipChain = LLHoc[ThisBin];
    int ipLast = ipChain;
    if (ipChain == -1){
      LLHoc[ThisBin] = ip;
    } else {
      ipLast = ipChain;
      while(ipChain != -1){
	if(particles[ipChain]->Dist8() > particles[ip]->Dist8())
	  break;
	ipLast = ipChain;
	ipChain = LL[ipChain];
      }
      LL[ip] = ipChain;
      if (ipLast == LLHoc[ThisBin]){
	LLHoc[ThisBin] = ip;
      } else {
	LL[ipLast] = ip;
      }
    }
  }
}
#endif

void MakeLL(int LLBins, int * &NInBin, int * &BinStart, vector <Particle *> &particles)
{
  //Determine the numbeer of bins
  LLBins = GalaxyZBin(ZREDMAX)+1;
  cout<<"    The list will have "<<LLBins<<" bins."<<endl;
  NInBin = new int[LLBins];
  BinStart = new int[LLBins];
  for(int i=0;i<LLBins;i++){
    NInBin[i] = 0;
    BinStart[i] = 0;
  }

  //Sort the particles in the z-bins in order of density
  int FirstParticle, LastParticle;
  float z0, z1;
  FirstParticle = 0;
  z0 = 0.0;
  for(int ib=0;ib<LLBins;ib++){
    cout<<"Sorting bin "<<ib<<endl;
    BinStart[ib] = FirstParticle;
    z1 = z0+(1+z0)*zTol;
    ip = FirstParticle;
    while(particles[ip]->Zred() < z1){
      ip++;
    }
    LastParticle = ip - 1;
    NInBin[ib] = LastParticle - FirstParticle + 1;
    sort(particles[FirstParticle], particles[LastParticle], DLessPart);
    z0 = z1;
    FirstParticle = LastParticle;
  }
}

float SelectGalaxyZ()
{
  //selects a random radius for galaxy s.t. prob /propto r^3
  float rn = drand48();
  rn = rn*rn*rn;
  rn *= cosmo.RofZ(RAMAX);
  float z = cosmo.ZofR(rn);
  return z;
}

void AssignGalaxies(vector <Particle *> &particles, vector <Galaxy *> &galaxies, ChainEl chel)
{
  // Sort galaxies on desired local mass density
  sort(galaxies.begin(),galaxies.end(),DLessGal);
  //Make Linked List of particles based on density
  int * LL;
  LL = new int[particles.size()];
  int LLBins = 0;
  int * LLHoc;
  cout<<"  Making the particle Linked List."<<endl;
  //MakeLL(LLBins, LLHoc, LL, particles);
  int * NInBin;
  int * BinStart;
  MakeLL(LLBins, NInBin, BinStart, particles);
  int LLpi[LLBins];//how far in the density linked list have we progressed at each z bin
  for(int pi=0;pi<LLBins;pi++)
    //LLpi[pi] = LLHoc[pi];
    LLpi[pi] = 0;
  int NCentralsSkipped = 0;//Number of centrals we find, i.e. ones that were assigned earlier
  int NUnderdense = 0;//Number of galaxies for which we couldn't find a dense enough halo
  float MaxRhoError = 0.; //The maximum difference in Dist8
  float Percent = 0.0;

  cout<<"  Looping through the Galaxies."<<endl;
  for(int gi=0;gi<galaxies.size();gi++){
    if(float(gi)/float(galaxies.size())> Percent){
      cout<<"    "<<100*Percent<<"% done"<<endl;
      Percent += 0.1;
    }
    //Was this galaxy already assigned as a central?
    if(galaxies[gi]->Central()){
      NCentralsSkipped++;
      continue;
    }
    //Select a redshift for the galaxy
    float zGal = SelectGalaxyZ();
    //Find the redshift bin of this galaxy
    int zBin = GalaxyZBin(zGal);
    float zCut = zTol*(1+zGal);
    //Select Galaxy Density
    float mag = galaxies[gi]->Mr();
    mag = mag + zGal*(zGal - 0.1);
    float ThisMStar = Mstar + Q*(zGal - 0.1);
    float lum = pow(10.,-0.4*(mag+ThisMStar));
    FiveTuple fTup(chel.cmean(),chel.fmean(lum),chel.fsig(lum),chel.ffrac(lum)); 
    galaxies[gi]->Dist8(fTup.LocalDens());

    //loop through the z bins and find out particle
    float dRho = 1e10; //How close are we to the desired density?
    int GalaxyParticle; //The particle we're keeping 
    for(int iz=zBin-1;iz<=zBin+1;iz++){
      if(iz < 0 || iz >= LLBins)
	continue;
      //int pi = LLpi[iz];
      int pi = LLpi[iz] + BinStart[iz];
      int piLast = 0;
      float dRhoLast = 0.0;
      while(pi != -1){
	if(abs(particles[pi]->Zred() - zGal) < zCut){
	  if(particles[pi]->Dist8() > galaxies[gi]->Dist8())
	    break;
	  dRhoLast = abs(particles[pi]->Dist8()-galaxies[gi]->Dist8());
	  piLast = pi;
	}
	pi = LL[pi];
      }
      //Make sure we found a particle in this bin
      if(pi == -1){
	NUnderdense++;
      } else {
	//find out of pi or piLast is closer in density
	float ThisdRho = abs(particles[pi]->Dist8()-galaxies[gi]->Dist8());
	if (dRhoLast < dRho){
	  ThisdRho = dRhoLast;
	  pi = piLast;
	}
	//Do we want to keep our best-fit galaxy from this bin?
	if(ThisdRho < dRho){
	  dRho = ThisdRho;
	  GalaxyParticle = pi;
	}
	//Update our LL counter
	LLpi[iz] = pi;
      }
    }
    //Now we can actually assign the galaxy
    galaxies[gi]->P(particles[GalaxyParticle]);
    if (dRho > MaxRhoError)
      MaxRhoError = dRho;
  }


  cout<<"First loop through galaxiees"<<endl;
  cout<<"Maximum RhoError: "<<MaxRhoError<<endl;
  cout<<"Skipped "<<NCentralsSkipped<<" previously assigned centrals."<<endl;
  cout<<NUnderdense<<" galaxies had particles that were too underdense"<<endl;

  //randomly assign any galaxies that weren't assigned above
  for(int gi=0;gi<galaxies.size();gi++){
    Particle * myp = galaxies[gi]->P();
    if((!myp)&&(galaxies[gi])){
      int pind= randint(particles.size()-1);
      galaxies[gi]->P(particles[pind]);
    }
  }

  //cout<<galaxies.size()<<"...";
  for_each(galaxies.begin(),galaxies.end(),DeleteAndNullifyBadgal);
  galaxies.erase(remove(galaxies.begin(),galaxies.end(),static_cast<Galaxy*>(0)), galaxies.end());
  //cout<<galaxies.size()<<endl;
  //cout<<"done."<<endl;
  cout<<"Finished assigning galaxies."<<endl;
}

vector <Particle *> AssignDimGalaxies(vector <Particle *> &particles, vector <Galaxy *> &galaxies)
{
    float pzmin = ZREDMAX;
    float pzmax = ZREDMIN;
    for(int pi=0;pi<particles.size();pi++){
      if(particles[pi]->Zred()<pzmin) 
	pzmin = particles[pi]->Zred();
      if(particles[pi]->Zred()>pzmax) 
	pzmax = particles[pi]->Zred();
    }
    //    cout<<"ag pzmin:"<<pzmin<<" pzmax:"<<pzmax<<" "<<particles.size()<<endl;


  //MSG("Assigning gals");
  //#define NOBIAS
  //MSG("[assigngals] Sorting gals");
  // Sort galaxies on desired local mass density
  //  vector <Particle *> particles;
  sort(galaxies.begin(),galaxies.end(),DLessGal);
  sort(particles.begin(),particles.end(),DLessP);
  int ng = galaxies.size();
  int np = particles.size();
  for(int gi=0;gi<galaxies.size();gi++){
    float best = gi*1.0/ng*np;
    int mini = (int) best-50;
    if (mini <=0) mini = 0;
    int maxi = (int) best+50;
    if(maxi<100) maxi = 100;
    if (maxi >=np) maxi = np-1;
    int ri;
    //    cout<<mini<<" "<<maxi<<" "<<gi<<" "<<np<<endl;
    ri = randint(mini, maxi);
    //   float g8 = galaxies[gi]->Dist8();
    float d8 = particles[ri]->Dist8();

    Point p = particles[ri]->BoxPosition();
    vector <float> deltax(3);
    bool outside = true;
    int j=0;
    while((outside)&&(j<1000)){
      for(int i = 0; i<3;i++){
	double r1 = (drand48()*d8*2-d8)/sim.LengthUnit(); 
	deltax[i] = r1;
      }
      //p.Print();
      p.Reset(p.X()+deltax[0], p.Y()+deltax[1], p.Z()+deltax[2]);
      float _ra = angle_const*atan2(p.Y(),p.X()); 
      if(_ra<0) _ra = _ra+360.;
      float _dec = angle_const*asin(p.Z()/p.R());
      if((_ra>RAMIN)&&(_ra<RAMAX)&&(_dec>DECMIN)&&(_dec<DECMAX))
	outside = false;
      //else{
      //	if(j==99)
      //cout<<j<<" "<<particles[ri]->Ra()<<" "<<particles[ri]->Dec()<<" ..."
      //    <<RAMAX<<" "<<DECMAX<<" "
      //    <<_ra<<" "<<_dec<<" ";p.Print();
      //      }
      j++;
    }

    Point v = particles[ri]->Velocity();
    Particle * particle = new Particle(p, v, d8);  
    //    if(j>99) {p.Print(); particles[ri]->PosPrint();}
    galaxies[gi]->P(particle);
    //    if(galaxies[gi]->Z()<pzmin) cout<<galaxies[gi]->Z()<<" "<<particles[ri]->Zred()<<" "<<particle->Zred()<<endl;
    //if(galaxies[gi]->Z()>pzmax) cout<<galaxies[gi]->Z()<<" "<<particles[ri]->Zred()<<" "<<particle->Zred()<<endl;
				    
    particles.push_back(particle);
  }
  return particles;
  
}


vector <Particle *> AssignGalaxiesRandomly(vector <Particle *> &particles, vector <Galaxy *> &galaxies)
{
  // Choose the galaxies at random from the particles
  vector <Particle *> fakeparticles;
  for(int gi=0;gi<galaxies.size();gi++){
    int ri;
    ri = randint(particles.size())-1;
    float d8 = particles[ri]->Dist8();
    Point p = particles[ri]->BoxPosition();
    vector <float> deltax(3);
    bool outside = true;
    for(int j=0;j<10;j++){
      for(int i = 0; i<3;i++){
	double r1 = (drand48()*d8*2-d8)/sim.LengthUnit(); 
	deltax[i] = r1;
      }
      //p.Print();
      p.Reset(p.X()+deltax[0], p.Y()+deltax[1], p.Z()+deltax[2]);
      if((p.X()>0)&&(p.Y()>0)&&(p.Z()>0)&&(p.X()<XMAX)&&(p.Y()<YMAX)&&(p.Z()<ZMAX)) {outside = false; break;}
    }
    //cout<<"delta"<<d8<<" "<<d8/sim.LengthUnit()<<" "<<deltax[0]<<" "<<deltax[1]<<" "<<deltax[2]<<endl;
    Point v = particles[ri]->Velocity();
    Particle * particle = new Particle(p, v, d8);  
    galaxies[gi]->P(particle);
    fakeparticles.push_back(particle);
  }
  return fakeparticles;
}

void AssignGalaxiesUniquely(vector <Particle *> &particles, vector <Galaxy *> &galaxies)
{
  // Choose the galaxies at random from the particles
  for(int gi=0;gi<galaxies.size();gi++){
    int ri;
    bool gotit = false;
    while(!gotit){
      ri = randint(particles.size())-1;
      if(!particles[ri]->IsGal()){
	gotit = true;
	//particles[ri]->MakeGal(galaxies[gi]->Gid());
	galaxies[gi]->P(particles[ri]);
      }
    }
  } 
}

//
// Note: There was an AssignPTGalaxies in here at one point... deleted in v 1.8
// Also deleted the old AssignGalaxies(particles, gals), but this should have been
// the same. (was used in pdfxi.cc
//

