#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include "biniostream.h"
#include "particle.h"
#include "halo.h"
#include "point.h"
#include "singleton.h"
#include <vector>
#include <list>

extern Simulation sim;

void ReadWHids(vector <Halo *> halos, vector <Particle *> particles){
  for(int pi=0;pi<10;pi++){
    particles[pi]->PosPrint();
    int tmp_id = particles[pi]->Hid();
    halos[tmp_id]->Print();
    halos[tmp_id-1]->Print();
    halos[tmp_id+1]->Print();
  }
}

//
// reads all halos from file and returns vector of halo pointers
//
vector <Halo*> ReadWHalos(void){
  vector <Halo*> halos;
  halos.reserve(1000000);
  string filename;
  //if(box == W1a)
  // filename = Dir()+"hvd1.m12";
  //else if(box == W1b)
  //filename = Dir()+"hvd1e100.m12";
  //else if(box == W2)
  //filename = Dir()+"hvd2.m12";
  if(box == W384a)
    filename = Dir()+"eb14.r8.m12";
  else if(box == W384b)
    filename = Dir()+"eb15.r8.m12";
  else if(box == W384c)
    filename = Dir()+"eb16.r8.m12";
  else{
    cerr<<"Inconsistent box_type"<<simulation<<" "<<box<<endl;    exit(1);
  }
    
  ifstream file(filename.c_str());
  if (file.fail()) {
    cerr<<"error: cannot open "<<filename<<endl;    exit(1);
  }
  else cerr<<"reading file "<<filename<<endl;
  //read in first two lines of strings.
  //for(int i=0; i<13;i++){
  //string tmps;
  //file>>tmps;
  //}
  int hid = 0;
  Point pos(0,0,0);
  Point vel(0,0,0);
  //the first halo is a false one
  Halo * halo = new Halo(pos, vel, 0, 0, hid, 0, 0, 0);
  halos.push_back(halo); hid++;

  //float junk;
  cout<<halos.size()<<endl;
  while(file){
    float m15,zred,sigma,xx,yy,zz,vx,vy,vz,rdelta;
    int pid; 
    int ip=0;
    file>>m15>>sigma>>xx>>yy>>zz>>vx>>vy>>vz>>pid>>rdelta;
      //	>>junk>>junk>>junk>>junk>>junk>>junk>>junk>>junk;
    assert(m15<3.5);
    assert(m15>sim.ParticleMass()*10./1e15);
    //    PRNTV(m15);
    //   cout<<"[read_warren] mh check:"
    //<<hid<<" "<<m15<<" "<<xx<<endl;
    Point pos(xx,yy,zz); 
    Point vel(vx,vy,vz);   
        
    zred = 0.0; //fix this
    //siglos = 0; I don't know why this was set to zero.


    
    Halo * halo = new Halo(pos, vel, m15*1e15, zred, ip, hid, sigma, rdelta);
    bool keep = true;
    //if(hid>0){if((xx==halos[hid-1]->X())&&
    //		 (yy==halos[hid-1]->Y())&&
    //	 (zz==halos[hid-1]->Z())) keep = false;}
    //    if(halo->M()<2.99e13) keep = false;
    //    if(halo->M()<=ParticleMass()*1e12*(nhpmin-1)) keep = false;
    //    if (halo->M()>1e15) cout<<halo->M()<<endl; else keep = false;
    if(keep){
      halos.push_back(halo);
      //cout<<hid<<" "<<halos.size()<<" "<<halo->Id()<<" "<<halos[halos.size()-1]->Id()<<" "<<halos[hid]->Id()<<endl;
      hid++;
    }
  }
  file.close();
  halos.pop_back();
  cerr<<"[read_warren] Read "<<halos.size()-1<<" halos."<<endl;
  cout<<halos[0]->X()<<" "<<halos[halos.size()-1]->X()<<" "<<halos[halos.size()-2]->X()<<endl;
  return halos;
}

//
// reads particles from file and returns vector of particle pointers
//
//vector <Particle *> ReadWarren(int &nread){
vector <Particle *> ReadWarren(int &nread, vector <Halo *> halos){
  vector <Particle*> particles;
  string pfilename, hfilename;
  /*  if(box == W1a){
    pfilename = Dir()+"lcdm_hvd1.1338.pmhdf";
    hfilename = Dir()+"hvd1.r8.txt";
  }
  else if(box == W1b){
    pfilename = Dir()+"lcdm_hvd1e100.0692.pmhdf";
    hfilename = Dir()+"hvd1e100.r8.txt";
  }
  else if(box == W2){
    pfilename = Dir()+"lcdm_hvd2.1338.pmhdf";
    hfilename = Dir()+"hvd2.r64.txt";
    }*/
  if(box == W384a){
    pfilename = Dir()+"lcdm_eb14.pmhdf";
    hfilename = Dir()+"eb14.r27.txt";
  }
  else if(box == W384b){
    pfilename = Dir()+"lcdm_eb15.pmhdf";
    hfilename = Dir()+"eb15.r27.txt";
  }
  else if(box == W384c){
    pfilename = Dir()+"lcdm_eb16.pmhdf";
    hfilename = Dir()+"eb16.r27.txt";
  }
  else{
    cout<<"Inconsistent box/sim type here"<<simulation<<" "<<box<<endl;exit(1);
  }

  binifstream pfile(pfilename.c_str());
  ifstream hfile(hfilename.c_str());
  //ifstream hidfile(fstrhid);

  //confirm that files open
  if (pfile.fail()) {
    cerr<<"error: cannot open particle file '" <<pfilename<<"'"<<endl;exit(1);
  }
  else cerr<<"reading file '"<<pfilename<<"'"<<endl;

  if (hfile.fail()) {
    cerr<<"error: cannot open halo file '" <<hfilename<<"'"<<endl;
    exit(1);
  }
  else cerr<<"reading file '"<<hfilename<<"'"<<endl;

#ifdef BIGENDIAN
  pfile.bigendian();
#endif

#ifdef CATALOG 
  string hidfilename = Dir()+"ahid."+MakeString(box,1)+".1.dat";
  ifstream hidfile(hidfilename.c_str()); 
  if (hidfile.fail()) {
    cerr<<"error: cannot open halo id file '" <<hidfilename<<"'"<<endl;
    exit(1);
  }
  else cerr<<"reading file '"<<hidfilename<<"'"<<endl;
#endif

  float tmp_x, tmp_y, tmp_z, dist8;
  float tmp_vx, tmp_vy, tmp_vz, tmp, tmp2;
  
  unsigned int np = 0;
  unsigned int nr = 0;


  particles.reserve(56623104);
#ifdef HALOASSIGN  
    while((pfile)&&(hfile)){
#endif
#ifdef CATALOG 
  int hid;
  while((pfile)&&(hfile)&&(hidfile)){
#endif
   //** Read 3 floats from particle file
    //** Convert to box positions as they are passed to constructor
    pfile>>tmp_x>>tmp_vx
	 >>tmp_y>>tmp_vy
	 >>tmp_z>>tmp_vz
	 >>tmp>>tmp2;
    //    cout<<tmp_x<<" "<<tmp_y<<" "
    //<<tmp_z<<"\t";
    //	//<<tmp_vx<<" "
    // <<tmp_vz<<" "<<tmp_vz<<"...  "
    // <<tmp<<" "<<tmp2<<" ";

    const double h_100 = 0.7;
    Point xx((tmp_x*h_100)/sim.Boxsize(), (tmp_y*h_100)/sim.Boxsize(), (tmp_z*h_100)/sim.Boxsize());
    //Point xx((tmp_x*h_100), (tmp_y*h_100), (tmp_z*h_100));
    Point vv(tmp_vx, tmp_vy, tmp_vz);

    //    if((xx.X()<0.0)||(xx.X()>1)) xx.Print();
    //if((xx.Y()<0.0)||(xx.Y()>1)) xx.Print();
    //if((xx.Z()<0.0)||(xx.Z()>1)) xx.Print();


    //** Read 1 short int from density file

    hfile>>dist8;//tmp_dens;
      
     
    Particle * particle = new Particle(xx, vv, dist8);
#ifdef CATALOG 
    hidfile>>hid;   //fix this
    if((hid>-1)&&(particle->Distance(halos[hid]->Position())>1.5*halos[hid]->R200()))
      cout<<"[read_warren]"<<hid<<" "
	  <<particle->Distance(halos[hid]->Position())<<" "
	  <<halos[hid]->R200()<<" "
	  <<particle->X()<<" "
	  <<halos[hid]->X()<<" "
	  <<endl;
    if(particle->Save()){
      particle->Hid(hid);
      particles.push_back(particle);
      np++;
      if(particle->Zred()>zmax->GetVal())
	zmax->SetVal(particle->Zred());    
    }
    else delete particle;
    if(nr>56623100){
      //what is this?
      cout<<"[read_warren: high particle number] "
	  <<nr<<" "<<hid<<"\t";
      cout<<dist8<<"\t";
      xx.Print();
    }
#endif
#ifdef HALOASSIGN
    particles.push_back(particle);
    np++;
#endif
    nr++;
    if(nr%5000000==0) PRNT("read_warren",nr);
    if(nr%5000000==0) PRNT(" . . . ",np);

  }
  if(pfile){
    cerr<<"pfile remains"<<endl;
    float tmp_float;
    int next = 0;
    while(pfile){
      pfile>>tmp_float>>tmp_float>>tmp_float>>tmp_float>>tmp_float>>tmp_float>>tmp_float>>tmp_float;
      next++;
    }
    cerr<<"n_ext="<<next<<endl;
  }
  if(hfile){
    cerr<<"hfile remains"<<endl;
    float tmp_float;
    int next = 0;
    while(hfile){
      hfile>>tmp_float;
      next++;
    }
    cerr<<"n_ext="<<next<<endl;
  }
#ifdef CATALOG 
  if(hidfile){
    cerr<<"hidfile remains"<<endl;
    int tmp_int;
    int next = 0;
    while(hidfile){
      hidfile>>tmp_int;
      next++;
    }
    cerr<<"n_ext="<<next<<endl;
  }
#endif
  //  if((pfile)||(hfile)||(hidfile)){
  // cerr<<pfile<<" "<<hfile<<" "<<hidfile<<endl;
  // cerr<<"[read_cube] Problem: files are different lengths"<<endl;
  // exit(1);
  //}
  cerr<<"[read_warren] "<<nr<<" particles read in this cube"<<"...";
  cerr<<"[read_warren] "<<np<<" particles saved."<<endl;
  assert(particles.size()>0);
  nread = nr;
  return particles;
}
