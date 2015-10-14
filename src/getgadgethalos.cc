#include <vector>
#include <iostream>
#include "particle.h"
#include "halo.h"
#include "biniostream.h"


//-----------> WARNING!  Currently this only works if you have a single gadget file<----------\\

vector <Particle *> ReadParticles(int &nread);

void HaloLinkedList(vector <Halo *> halos, int llbins, vector <int> &ll, vector <int> &llhoc){
  float scalefac;
  int ix, iy, iz, ihoc;
  scalefac = float(llbins)/float(sim.LengthUnit());
  cout<<"scalefac = "<<scalefac<<endl;
  for(int ih=0;ih<halos.size();ih++){
    //cout<<ih<<endl;
    ix = int(halos[ih]->X()*scalefac);
    iy = int(halos[ih]->Y()*scalefac);
    iz = int(halos[ih]->Z()*scalefac);
    ihoc = ix + llbins*iy + llbins*llbins*iz;
    ll[ih] = llhoc[ihoc];
    llhoc[ihoc] = ih;
    //    if (ih == 0 || ih == 1 || ih == 2)
    //      cout<<"Putting a halo in to chain "<<ix<<" "<<iy<<" "<<iz<<endl;
  }
}

void GetGadgetHalos(vector <Halo *> halos){
  int nread;
  vector <Particle*> particles = ReadParticles(nread);
  string fstrhid = Dir()+"analysis/group_info/dataoutput/ahid."+simlabel;
  ofstream pfile (fstrhid.c_str());
  vector <int> closest_halo(particles.size());
  for(int ih=0;ih<particles.size();ih++)
    closest_halo[ih]=-1;
  cout<<"Writing file."<<endl;
  float next_percent = 0.0;

  int llbins = 128;
  vector <int> ll(halos.size());
  vector <int> llhoc(llbins*llbins*llbins);
  for(int ill=0;ill<halos.size();ill++)
    ll[ill] = -1;
  for(int ix=0;ix<llbins*llbins*llbins;ix++)
    llhoc[ix] = -1;

  cout<<"Making the linked list..."<<endl;
  HaloLinkedList(halos, llbins, ll, llhoc);

  cout<<"Traversing the linked list...."<<endl;
  float scalefac = float(llbins)/float(sim.LengthUnit());
  for(int ip=0;ip<particles.size();ip++){
    if (float(ip)/float(particles.size()) >= next_percent){
      cout<<next_percent*100.<<"% done..."<<endl;
      next_percent += 0.01;
    }
    int ix = particles[ip]->X()*scalefac;
    int iy = particles[ip]->Y()*scalefac;
    int iz = particles[ip]->Z()*scalefac;
    //    if(ip == 0){
    //      cout<<"Particle Position = "<<particles[ip]->X()<<" "<<particles[ip]->Y()<<" "<<particles[ip]->Z()<<endl;
    //      cout<<"Particle cell: "<<ix<<" "<<iy<<" "<<iz<<endl;
    //    }
    double min_dist = 100.0;
    double AssignedM200 = 0.;
    int AssignedHalo = 0;
    closest_halo[ip] = -1;
    for(int ixx=ix-1;ixx<=ix+1;ixx++){
      if(ixx < 0 || ixx >= llbins) continue;
      for(int iyy=iy-1;iyy<=iy+1;iyy++){
	if(iyy < 0 || iyy >= llbins) continue;
	for(int izz=iz-1;izz<=iz+1;izz++){
	  if(izz < 0 || izz >= llbins) continue;
	  int lhoc = ixx + llbins*iyy + llbins*llbins*izz;
	  int lnext = llhoc[lhoc];
	  while(lnext >= 0){
	    //	    if (ip == 0){
	    //	      cout<<"Halo id: "<<halos[lnext]->Id()<<endl;
	    //	      cout<<"Halo Position: "<<halos[lnext]->X()<<" "<<halos[lnext]->Y()<<" "<<halos[lnext]->Z()<<endl;
	    //	    }
	    float dx = particles[ip]->X() - halos[lnext]->X();
	    float dy = particles[ip]->Y() - halos[lnext]->Y();
	    float dz = particles[ip]->Z() - halos[lnext]->Z();
	    float dist=sqrt(dx*dx + dy*dy + dz*dz);
	    if(dist < halos[lnext]->R200()){
	      if(halos[lnext]->M() > AssignedM200){
		min_dist = dist;
		closest_halo[ip] = halos[lnext]->Id();
		AssignedM200 = halos[lnext]->M();
		AssignedHalo = 1;
	      }
	    }
	    else if(dist < min_dist && AssignedHalo == 0){
	      min_dist = dist;
	      closest_halo[ip]=halos[lnext]->Id();
	    }
	    lnext = ll[lnext];
	  }
	}
      }
    }
	    
    /*
    if(ip < 10){
      cout<<"Closest halo id = "<<closest_halo[ip]<<", "<<min_dist<<endl;
      cout<<"   Particle Position: "<<particles[ip]->X()<<" "<<particles[ip]->Y()<<" "<<particles[ip]->Z()<<endl;
      cout<<"   Halo Position: "<<halos[closest_halo[ip]]->X()<<" "<<halos[closest_halo[ip]]->Y()<<" "<<halos[closest_halo[ip]]->Z()<<endl;
    }
    */
    pfile<<closest_halo[ip]<<endl;
  }
}

