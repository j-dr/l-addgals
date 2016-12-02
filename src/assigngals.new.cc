#include <vector>
#include <cassert> 
#include <iostream>
#include "galaxy.h"
//#include "choose.h"
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
bool MLessGal(Galaxy * a, Galaxy * b)
{
  return a->Mr() < b->Mr();
}

// compare galaxy pointers by local mass density
bool DLessP(Particle * a, Particle * b)
{
  return a->Dist8() < b->Dist8(); 
}

struct CandidatePart{
  int pi;
  float zred;
  float d8;
};

bool ClosestZ(CandidatePart a, CandidatePart b)
{
  return a.zred < b.zred;
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

void MakeLL(int &LLBins, int * &NInBin, int * &BinStart, vector <Particle *> &particles)
{
  //Determine the numbeer of bins
  LLBins = GalaxyZBin(ZREDMAX)+1;
  cout<<"    The list will have "<<LLBins<<" z-bins."<<endl;
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

  //first we need to determine the starting bin
  int startbin = 0;
  float old_z0;
  old_z0 = 0.0;
  while(z0 < particles[0]->Zred()){
    old_z0 = z0;
    z1 += z0+(1+z0)*zTol;
    z0 = z1;
    startbin++;
  }
  if(startbin > 0) 
    startbin--;
  z0 = old_z0;

  for(int ib=startbin;ib<LLBins;ib++){
    BinStart[ib] = FirstParticle;
    z1 = z0+(1+z0)*zTol;
    //cout<<"Bin "<<ib<<" covers the range z = "<<z0<<"-"<<z1<<endl;
    int ip = FirstParticle;
    while(ip < particles.size() && particles[ip]->Zred() < z1){
      //if(particles[ip]->IsGal()) //All particles added for centrals are at the end
	//break;                   //Note: This should do nothing.  
      ip++;
      if(ip >= particles.size())
	break;
    }
    LastParticle = ip;
    NInBin[ib] = LastParticle - FirstParticle;
    sort(&particles[FirstParticle], &particles[LastParticle], DLessP);
    if (FirstParticle >= 0 && FirstParticle <= particles.size()-3){
      cout<<"  Bin "<<ib<<" from z = "<<z0<<" to "<<z1<<":"<<endl;
      cout<<"    First Few z's:  "<<particles[FirstParticle]->Zred()<<" "<<particles[FirstParticle+1]->Zred()<<" "<<particles[FirstParticle+2]->Zred()<<endl;
      cout<<"    First Few d8's: "<<particles[FirstParticle]->Dist8()<<" "<<particles[FirstParticle+1]->Dist8()<<" "<<particles[FirstParticle+2]->Dist8()<<endl;
    }
    z0 = z1;
    FirstParticle = LastParticle;
  }
  int NPartInList = 0;
  int MaxBin = 0;
  int MinBin = 10000000;
  float NAvg = 0.;
  float NSig = 0.;
  for(int i=0;i<LLBins;i++){
    NPartInList += NInBin[i];
    NAvg += float(NInBin[i]);
    NSig += float(NInBin[i]*NInBin[i]);
    if (NInBin[i] > MaxBin)
      MaxBin = NInBin[i];
    if (NInBin[i] < MinBin)
      MinBin = NInBin[i];
  }
  NAvg /= float(LLBins);
  NSig /= float(LLBins);
  NSig = sqrt(NSig - NAvg*NAvg);
  cout<<"Made the list with "<<NPartInList<<" particles."<<endl;
  cout<<"  List Statistics:  Avg = "<<NAvg<<" Sigma = "<<NSig<<" min/max = "<<MinBin<<"/"<<MaxBin<<" Nbins = "<<LLBins<<endl;
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

