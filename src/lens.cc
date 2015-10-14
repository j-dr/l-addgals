#include "hv.h"
#include "cosmo.h"
#include "functions.h"
#include "cluster.h"
#include "box.h"

static const float zcut = 0.28;
//enum lens_center_type{HALOS, CLUSTERS};
//const lens_center_type lens_center = HALOS;


#ifdef HVL
Simulation sim(box, use_cells);
#endif
#ifdef LANL
Simulation sim(box);
#endif

Cosmology cosmo = sim.SimCosmology();


//#define CLUSTERS
#define HALOS

int main(void){
  string label = simlabel+"_"+flabel;
  cout<<sim.ParticleMass()<<endl;
  cout<<"dir:"<<Dir()<<endl;

#ifdef CLUSTERS
    label = label+"clust";
#endif
  cout<<"writing to label"<<label<<endl;

  //** Read in z as a function of distance so that redshifts 
  //** can be determined without integrating

  cosmo.ReadZFile();
  ReadKernelFile();
  //** Read in particles from cubes
  int nread =0;
  vector <Halo*> halos1 = ReadHalos();
  cout<<"Read "<<halos1.size()<<" halos."<<endl;

  vector <Particle*> particles = ReadParticles(nread, halos1);
  cout<<"Done."<< particles.size()<<endl;
  assert(particles.size()>0);


  string outfh = "/data1/risa/lens/halos"+label+".dat";
  string outfn3D = "/data1/risa/lens/lens_3Dmass"+label+".dat";
  string outfn2D = "/data1/risa/lens/lens_2Dmass"+label+".dat";
  string outfn2Datr = "/data1/risa/lens/lens_2Dmassatr"+label+".dat";
  string outfn3Drho = "/data1/risa/lens/lens_3Ddensity"+label+".dat";
  string outfnDS = "/data1/risa/lens/lens_DeltaSigma"+label+".dat";
  string outfnDSZC = "/data1/risa/lens/lens_DeltaSigmaZcut"+label+".dat";
  string outfnDSNW = "/data1/risa/lens/lens_DeltaSigmaNoWeight"+label+".dat";

  ofstream outh(outfh.c_str());
  ofstream out3D(outfn3D.c_str());
  ofstream out3Drho(outfn3Drho.c_str());
  ofstream out2D(outfn2D.c_str());
  ofstream out2Datr(outfn2Datr.c_str());
  ofstream outDS(outfnDS.c_str());
  ofstream outDSZC(outfnDSZC.c_str());
  ofstream outDSNW(outfnDSNW.c_str());
  vector <float> rmin, rmid;
  float binsize = 0.1;  //log bins
  float maxdist = 10.; //max search radius in Mpc/h
  for(double rr = -1.025;rr<=log10(maxdist);rr+=binsize){
    double ten = 10.;
    rmin.push_back(pow(ten,rr));
    rmid.push_back(pow(ten,(rr+0.5*binsize)));
  }
  int nb = rmin.size();

  ofstream outtest("lens.out");  
  ofstream fractest("lens.out");  
  int hkeep = 0;
  vector <float> sigltr(nb),sigatr(nb),sigltr_zcut(nb),sigatr_zcut(nb);
  vector <float> unwltr(nb),unwatr(nb),three_d(nb), three_datr(nb);

  vector <float> padist; //projected angular distance from the halo in comoving Mpc/h
  vector <float> ppdist; //3-D distance from the halo in comoving Mpc/h
  vector <float> pzz;

#ifdef CLUSTERS
    vector <Cluster *> halos = ReadClusters();
    cout<<"Read "<<halos.size()<<" clusters."<<endl;
#endif
#ifdef HALOS
      vector <Halo *> halos = halos1;
#endif

      cout<<"np:"<<particles.size()<<" nh:"<<halos.size()<<endl;

  for(int hi=0;hi<halos.size();hi++){
#ifdef HALOS
    Halo * h = halos[hi];
    double haloz = h->ZredReal();
    double halor = h->R();
    if((haloz<zcut)
       &&(h->Dec()<DECMAX-1)
       &&(h->Dec()>DECMIN+1)
       &&(h->Ra()>RAMIN+1)
       &&(h->Ra()<RAMAX-1)
       //       &&(h->X()<187.5*use_cells)&&(h->Y()<187.5*use_cells)&&(h->Z()<187.5*use_cells)
       &&(h->X()>0)&&(h->Y()>0)&&(h->Z()>0)
       &&(h->M()>5e13)){
#endif
#ifdef CLUSTERS
    Cluster * h = halos[hi];
    double haloz = h->ZredReal();
    double halor = h->R();
    if((haloz<zcut)
       &&(h->Dec()<DECMAX-1)
       &&(h->Dec()>DECMIN+1)
       &&(h->Ra()>RAMIN+1)
       &&(h->Ra()<RAMAX-1)){
       //       &&(h->Ngals()>4)){
      cout<<"clust: "<<haloz<<" "<<h->Ra()<<" "<<h->Dec()<<" "<<h->Ngals()<<endl;
#endif
      hkeep++;
      double mindist = 100;
      Particle * minp = 0;
      int pp=0;
      for(int pi=0;pi<particles.size();pi++){
	Particle * p = particles[pi];
	// comoving distance on a sphere at z equals halo z
	double com_theta = h->Theta(p)*halor;
	double threeddist = p->Distance(h->Position());
	//need to convert the halo radius into comoving coordinates
	if(threeddist<h->R200()*(1+haloz)){
	  pp++;
	}
	if(threeddist<mindist){
	  mindist = threeddist;
	  minp = p;
	}
	if(hkeep==0)
	  fractest<<com_theta<<" "<<threeddist<<endl;
	if(threeddist<maxdist){
	  padist.push_back(com_theta);
	  ppdist.push_back(threeddist);
	  pzz.push_back(p->ZredReal());
	}
      }
      //      this is shrink to fit.
      cout<<"particles saved:"<<ppdist.size()<<" "<<ppdist.capacity()<<" ";
      vector <float>(ppdist).swap(ppdist);
      vector <float>(padist).swap(padist);
      vector <float>(pzz).swap(pzz);
      cout<<ppdist.size()<<" "<<ppdist.capacity()<<" ";
      
      if(!pp)
	cerr<<"halo out of bounds"<<pp<<" "<<mindist<<" "<<h->X()<<" "<<h->Y()<<" "<<h->Z()<<" "<<h->Ra()<<" "<<h->Dec()<<endl;
      else{
	cerr<<pp<<" "<<mindist<<" "<<h->X()<<" "<<h->Y()<<" "<<h->Z()<<" "
	    <<minp->X()<<" "<<minp->Y()<<" "<<minp->Z()<<endl;
	cerr<<pzz.size()<<" particles saved."<<endl;

#ifdef HALOS	
      outh<<hi<<" "<< h->M()<<" "<<h->R200()<<" "<<h->Zred()<<" "<<h->Sig()<<" "//endl;
	  <<pp<<" "<<mindist<<" "<<h->X()<<" "<<h->Y()<<" "<<h->Z()<<" "
	  <<minp->X()<<" "<<minp->Y()<<" "<<minp->Z()<<" "
	  <<h->Ra()<<" "<<h->Dec()<<" "<<endl;
#endif

#ifdef CLUSTERS	
      outh<<hi<<" "<< h->Ngals()<<" "<<h->Zred()<<" "
	  <<pp<<" "<<mindist<<" "<<h->X()<<" "<<h->Y()<<" "<<h->Z()<<" "
	  <<minp->X()<<" "<<minp->Y()<<" "<<minp->Z()<<" "
	  <<h->Ra()<<" "<<h->Dec()<<endl;
#endif
      for(int i=0;i<nb;i++){
	sigltr[i]=0;sigatr[i]=0;sigltr_zcut[i]=0;sigatr_zcut[i]=0;
	three_d[i]=0.;unwatr[i]=0;unwltr[i]=0; three_datr[i]=0.;
	
	for(int ai=0;ai<padist.size();ai++){
	  if((padist[ai]>=rmin[i])&&(padist[ai]<rmin[i+1])){
	    sigatr[i] += SigCritInv(pzz[ai])*(sim.ParticleMass()/1e12); //weight*particle mass in 1e12 Msun/h
	    unwatr[i] += sim.ParticleMass()/1e12; 	    //unweighted*particle mass in 1e12 Msun/h
	    if(pzz[ai]<zcut){
	      sigatr_zcut[i] += SigCritInv(pzz[ai])*(sim.ParticleMass()/1e12); //weight*particle mass in 1e12 Msun/h/pi
	    }
	  }
	  if(padist[ai]<rmid[i]){
	    sigltr[i] += SigCritInv(pzz[ai])*(sim.ParticleMass()/1e12); //weight*particle mass in 1e12 Msun/h
	    unwltr[i] += sim.ParticleMass()/1e12;
	    if(pzz[ai]<zcut){
	      sigltr_zcut[i] += SigCritInv(pzz[ai])*(sim.ParticleMass()/1e12); //weight*particle mass in 1e12 Msun/h
	    }
	    if(ppdist[ai]<rmid[i]){
	      //mass enclosed in a 3-D sphere.
	      three_d[i] += sim.ParticleMass();
	    }
	    if((ppdist[ai]>=rmin[i])&&(ppdist[ai]<rmin[i+1])){
	      //mass enclosed in a 3-D sphere.
	      three_datr[i] += sim.ParticleMass();
	    }
	    //cout<<i<<" "<<notcount<<" not counted for low z"<<endl;
	  }
	}//end of particle loop
      } //end bin loop
      
      //mass unit is 10^12 and size unit is 10^6
      //so this outputs in erin's units M_sun/pc^2
      for(int i=0;i<rmin.size()-1;i++){
	cout<<rmid[i]<<" ";
	sigltr[i] = sigltr[i]/SigCritInv(haloz);
	sigatr[i] = sigatr[i]/SigCritInv(haloz);
	sigltr_zcut[i] = sigltr_zcut[i]/SigCritInv(haloz);
	sigatr_zcut[i] = sigatr_zcut[i]/SigCritInv(haloz);
	outDS<<sigltr[i]/PI/sqr(rmid[i])-sigatr[i]/PI/(sqr(rmin[i+1])-sqr(rmin[i]))<<" ";
	outDSNW<<unwltr[i]/PI/sqr(rmid[i])-unwatr[i]/PI/(sqr(rmin[i+1])-sqr(rmin[i]))<<" ";
	//outDSNW<<unwatr[i]<<" ";
	outDSZC<<sigltr_zcut[i]/PI/sqr(rmid[i])-sigatr_zcut[i]/PI/(sqr(rmin[i+1])-sqr(rmin[i]))<<" ";
	out3D<<three_d[i]<<" ";
	out3Drho<<three_datr[i]/(4.*PI/3.*(rmin[i+1]*rmin[i+1]*rmin[i+1]-rmin[i]*rmin[i]*rmin[i]))<<" ";
	cout<<unwltr[i]<<" ";
	out2D<<unwltr[i]<<" ";
	out2Datr<<unwatr[i]<<" ";
      } 
      cout<<endl;
      outDS<<endl;
      outDSNW<<endl;
      outDSZC<<endl;
      out3D<<endl;
      out3Drho<<endl;
      out2D<<endl;
      out2Datr<<endl;
    ppdist.clear();
    padist.clear();
    pzz.clear();
    cout<<".."<<ppdist.size()<<" "<<ppdist.capacity()<<" ";
    //this is really clearing it.
    vector <float>().swap(ppdist);
    vector <float>().swap(padist);
    vector <float>().swap(pzz);
    cout<<ppdist.size()<<" "<<ppdist.capacity()<<endl;
      } //end if particle loop
    }
  } //end halo loop
} //end main

