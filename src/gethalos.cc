#include "box.h"
#include "hv.h"
//#include "stl_util.h" //for copy_if

#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <list>
#include <strstream>

#include "biniostream.h"
#include "point.h"
#include "particle.h"
#include "halo.h"
#include "singleton.h"
#include "functions.h"             


#ifdef HVL
Simulation sim(box, use_cells);
#endif
#ifdef LANL
Simulation sim(box);
#endif
#ifdef GADGET
Simulation sim(box);
#endif

Cosmology cosmo = sim.SimCosmology();
using namespace std;
void GetWHids(vector <Halo *> halos, vector <Particle *> particles);
void GetGadgetHalos(vector <Halo *> halos);

class CloseHalo{
public:
  CloseHalo(Point center, float within):_center(center),_within(within){}
  bool operator()(Halo* h) const;
  
private:
  Point _center;
  float _within;
};

inline bool CloseHalo::operator()(Halo* h) const {
  float distance =_center.Distance(h->Position());
  return (distance<_within);
}

Point Center(int i, int j, int k){
  float c;
  if((box==HVLIGHT)||(box==HVZ0)||(box==HVPOW)){
    c = 0.5;
  }
  else 
    c = 0.0;
  float x = ((i+c)/BOXES)*sim.LengthUnit();
  float y = ((j+c)/BOXES)*sim.LengthUnit();
  float z = ((k+c)/BOXES)*sim.LengthUnit();
  if(box==HVLIGHT){
    x -= 0.5*sim.LengthUnit();
    y -= 0.5*sim.LengthUnit();
    z -= 0.5*sim.LengthUnit();
  }
  cout<<"[center]"<<i<<" "<<j<<" "<<k<<" center: "<<x<<" "<<y<<" "<<z<<endl;
  Point p(x,y,z);
  return p;
}

int ReadCube(int i, int j, int k, vector <Particle *> &particles, 
	     int xstart, int ystart, int zstart);

int main(void){
  //read in all the halos
  vector <Halo* > halos = ReadHalos();
  PRNTV(halos.size());
  int isolated = 0;
  if(simulation==HV){
    //read in particles one cube at a time.
    int xstart = 4;
    int ystart = 7;
    int zstart = 7;
    int nxcubes = 1;
    int nycubes = 1;
    int nzcubes = 1;
    int xfinish = xstart+nxcubes-1;
    int yfinish = ystart+nycubes-1;
    int zfinish = zstart+nzcubes-1;
    
    float search_rad = 235.0;
    for(int i=xstart;i<=xfinish;i++) 
      for(int j=ystart;j<=yfinish;j++) 
	for(int k=zstart;k<=zfinish;k++){
	  vector <Halo*> shalos;
	  Point center=Center(i,j,k);
	  CloseHalo ch(center, search_rad);
	  // copy all the halos within search_rad of the center of the cube.
	  copy_if(halos.begin(),halos.end(),
		  back_inserter(shalos), ch);
	  
	  cout<<"looking at "<<shalos.size()<<" halos"<<endl;
	  if (shalos.size() > 0)
	    shalos[0]->Print();
	  center.Print();
	  for(int ig=0;ig<shalos.size();ig++){
	    if(shalos[ig]->Id() == 983)
	      cout<<" Including halo 983!"<<endl;
	    }
	  //** Read in particles from cubes
	  vector <Particle *> particles;
	  //** Reserve max number of particles to avoid memory difficulties
	  particles.reserve(MAXSIZE);
	  if(particles.size()>MAXSIZE){
	    cerr<<"Exceeded maximum number of particles"
		<<particles.size()<<" "<<MAXSIZE<<endl;
	    exit(1);
	  }
	  cout<<"Checking cube"<<i<<" "<<j<<" "<<k<<"..."<<endl;
	  ReadCube(i,j,k, particles, xstart, ystart, zstart);
	  cout<<particles.size()<<" particles saved."<<endl;
	  if (particles.size() > 0){
	    cout<<" First Particle Position Info: ";
	    particles[0]->Print();
	    particles[1]->Print();
	    particles[2]->Print();
	  }

	  if(particles.size()>0){
	  string fstrhid = Dir()+"rnn/ahid."+MakeString(i,2)+"."+MakeString(j,2)+"."+MakeString(k,2);	
	  ofstream pfile(fstrhid.c_str());
	    vector <int> closest_halo(particles.size());
	    cout<<"Writing file:"<<i<<" "<<j<<" "<<k<<endl;
	    
	    for(int pi=0;pi<particles.size();pi++){
	      vector <double> dist(shalos.size());
	      double min_dist = 100.0;
	      closest_halo[pi] = -1;
	      for(int hi=0;hi<shalos.size();hi++){
		//dist[hi] = particles[pi].Distance(halos[hi].Position());
		dist[hi] = particles[pi]->Distance(shalos[hi]->Position());
		if(dist[hi]<min_dist){
		  //cout<<dist[hi]<<endl;
		  min_dist = dist[hi];
		  closest_halo[pi] = shalos[hi]->Id();
		}
	      }
	      if(min_dist>=100) {//cout<<"*** no close halo: "; 
		isolated++;};
	      //	      if(min_dist>70) cout<<"big min dist:"<<min_dist<<" "<<i<<" "<<j<<" "<<k<<" "<<particles[pi]->X()<<" "
	      //		  <<particles[pi]->Y()<<" "
	      //		  <<particles[pi]->Z()<<" "
	      //	  <<endl;
	      //pfile<<closest_halo[pi]<<endl;
	      pfile<<closest_halo[pi]<<" "<<min_dist<<endl;
	    }
	  }
	}
    cout<<isolated<<"isolated particles found.  Done."<<endl;
  }
  else if(simulation==WAR){
    vector <Particle *> ReadWarren(int &nread, vector<Halo *> halos);
    int nread;
    vector <Particle *> particles = ReadWarren(nread, halos);
    PRNTV(nread);
    PRNTV(particles.size());
    MSG("getting halo ids");
    GetWHids(halos,particles);


#ifdef old_junk
    string fstrhid = Dir()+"ahid."+MakeString(box,1)+".2.dat";	
    ofstream pfile(fstrhid.c_str());
    
    vector <int> closest_halo(particles.size());
    
    for(int pi=0;pi<particles.size();pi++){
      int halo_id = -1;
      float hdist = 100;
      //      vector <double> dist(shalos.size());
      for(int hi=1;hi<halos.size();hi++){
	if(hi!=halos[hi]->Id()) cout<<hi<<" "<<halos[hi]->Id()<<endl;
	float dist = particles[pi]->Distance(halos[hi]->Position());
	if(dist<halos[hi]->R200()/1000.){
	  //assert(hi==halos[hi]->Id());
	  pfile<<hi<<endl;
	  break;
	}
	else if(dist<hdist){
	  hdist = dist;
	  halo_id = hi;
	}
      }
      pfile<<halo_id<<endl;
    }
#endif

  }
  else if(simulation==GADGET2){
    GetGadgetHalos(halos);
  }
    // if(particles[pi]->InHalo()){
    //cout<<pi<<"\t";      cout<<hid<<endl;//"\t";
    //particles[pi]->Print();
    //assert(hid<halos.size());
    //if(hid<halos.size()){
    //  if(hid>=0) 
    //    halos[hid]->Print();
    //  else{
    //    isolated++;
    //    cout<<"no halo\n";
    //  }
    //}
    //}
    //}
    //}
}
