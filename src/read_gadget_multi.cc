#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <list>
#include <strstream>
#include <iostream>
#include <sstream>
#include <string>

#include "biniostream.h"
#include "point.h"
#include "particle.h"
#include "halo.h"
#include "singleton.h"
#include "functions.h"
#include "constants.h"
#include "gadget.h"

static const int MMAXINT=32767;


//
// reads all halos from file and returns vector of halo pointers
//

//routine for objects in ROCKSTAR ascii format
vector <Halo*> ReadRockstarHalos(void){
  vector <Halo*> halos;
  string filename, rnnfilename;

  filename = halofile;
  rnnfilename = rnn_halofile;
  ifstream file(filename.c_str());
  if (file.fail()) {
    cerr<<"[ReadRockstarHalos error]: cannot open "<<filename<<endl;    exit(1);
  }
  ifstream rfile(rnnfilename.c_str());
  if (rfile.fail()) {
    cerr<<"[ReadRockstarHalos error]: cannot open "<<rnnfilename<<endl;    exit(1);
  }

  //#  id, did, m200c, vmax, vrms, r200c, rs, np, x, y, z, vx, vy, vz, pid

  //read in first line of strings.
  string tmps;

  int hid = 0;
  while(file){
    double mvir;
    float xMpc,yMpc,zMpc,vx,vy,vz,m200c, r200c, vmax, vrms, rs, zred, rdel;
    int ip, hid, pid, did;
    file>>hid>>did>>m200c>>vmax>>vrms>>r200c>>rs>>ip>>xMpc>>yMpc>>zMpc>>vx>>vy>>vz>>pid;
    rfile>>rdel;
    zred = cosmo.ZofR(sqrt(xMpc*xMpc + yMpc*yMpc + zMpc*zMpc));

    if (hid < 10){
      cout<<"Halo Position "<<hid<<" = "<<xMpc<<" "<<yMpc<<" "<<zMpc<<", host id = "<<pid<<endl;
    }

    //if(m200c < BCG_Mass_lim) continue;
    //if(pid >= 0) continue;
    Point vel(vx,vy,vz);
    float xx, yy, zz;

    xx = xMpc;
    yy = yMpc;
    zz = zMpc;

    Point pos(xx, yy, zz);
    Halo * halo = new Halo(pos, vel, m200c, zred, ip, hid, vrms, r200c);
    halo->Dist8(rdel);
    //halo->Particle(pid);
    halo->Host(pid);

#ifdef SNAPSHOT
    if(xx > 0 && yy > 0 && zz > 0 &&
       xx < sim.Boxsize()*BOXFR && yy < sim.Boxsize()*BOXFR &&
       zz < sim.Boxsize()*BOXFR)
      {
        halos.push_back(halo);
        hid++;
      }
#else
    halos.push_back(halo);
    hid++;
#endif

  }
  file.close();
  cout<<"[ReadRockstarHalos] Read "<<halos.size()<<" halos."<<endl;
  return halos;
}


//Routine for objects in Risa's original ascii format
vector <Halo*> ReadGadgetHalos(void){
  vector <Halo*> halos;

  string filename;
  filename = halofile;
  ifstream file(filename.c_str());
  if (file.fail()) {
    cerr<<"[ReadGadgetHalos error]: cannot open "<<filename<<endl;    
    exit(1);
  }
  cout<<"Opened file for reading halos: "<<filename<<endl;
    //Expected file format:
    //#  mvir    zred sigma ip    xMpc    yMpc    zMpc    vx    vy    vz siglos rdelta part_id

    int hid = 0;
    while(file){
      double mvir;
      float zred,sigma,xMpc,yMpc,zMpc,vx,vy,vz,siglos,rdelta, rdel;
      int ip, id, pid; 
      file>>id>>mvir>>sigma>>ip>>xMpc>>yMpc>>zMpc>>vx>>vy>>vz>>siglos>>rdelta>>rdel>>pid;
      zred = cosmo.ZofR(sqrt(xMpc*xMpc + yMpc*yMpc + zMpc*zMpc));
      
#ifdef DEBUG
      if (hid < 10){
	cout<<"Halo Position "<<hid<<" = "<<xMpc<<" "<<yMpc<<" "<<zMpc<<", particle at center = "<<pid<<endl;
      }
#endif

      if(mvir < BCG_Mass_lim) continue;
      Point vel(vx,vy,vz);   
      float xx, yy, zz;
    
      xx = xMpc;
      yy = yMpc;
      zz = zMpc;
    
      Point pos(xx, yy, zz);
      Halo * halo = new Halo(pos, vel, mvir, zred, ip, hid, siglos, rdelta);
      halo->Dist8(rdel);
      halo->Particle(pid);

#ifdef SNAPSHOT
      if(xx > 0 && yy > 0 && zz > 0 &&
	 xx < sim.Boxsize()*BOXFR && yy < sim.Boxsize()*BOXFR &&
	 zz < sim.Boxsize()*BOXFR)
	{
          halos.push_back(halo);
          hid++;
        }
#else
      halos.push_back(halo);
      hid++;
#endif     
 
    }
    file.close();

  cout<<"[read_cube] Read "<<halos.size()<<" halos."<<endl;
  return halos;
}

int ReadGadgetnFiles(void){
  int buf;
  string fstrp = simulationfile;

  //We need to determine which form of the gadget file exists, very inelegant
  const char *fcharp = fstrp.c_str();
  ifstream inFile(fcharp, ios::in);
  if(!(inFile.good()))
    {
      inFile.close();
      fstrp = fstrp+".0";
      const char *fcharp2 = fstrp.c_str();
      ifstream inFile2(fcharp2, ios::in);
      if(!(inFile2.good()))
	{
	  cout<<"Error!  The file `"<<fstrp<<"' does not exist w/ or w/o the trailing number"<<endl;
	  inFile2.close();
	  return 0;
	}
      inFile2.close();
    }
  inFile.close();

  binifstream pfile(fstrp.c_str());

  //confirm that files open
  if (pfile.fail()) {
      cerr<<"[ReadGadgetnFiles] error: cannot open file '" <<fstrp<<"'"<<endl;
      return 0;
  }

  gadget_header header;

  //read the gadget header -- most of the information we don't care about
  cout<<"Reading particle data from file "<<fstrp<<endl;
  pfile>>buf; //it's f77_unformatted
  if (buf != 256){
    cout<<"Incorrect header format for file "<<fstrp<<", buf = "<<buf<<endl;
    return 0;
    }
  pfile.read((char *)(&header), sizeof(struct gadget_header));

  pfile.close();
  return header.num_files;
}

int ReadGadget(int nfiles, vector <Particle *> &particles){
  gadget_header header;
  int buf;
  int nr = 0;
  int np = 0;

  //define out input files -- one for the particles, one for rnn
  string base_fstrp = simulationfile;
  string base_fstrr = rnnfile;

  float maxpos[3];
  float minpos[3];
  maxpos[0] = maxpos[1] = maxpos[2] = 0.;
  minpos[0] = minpos[1] = minpos[2] = 10000000.;

  nr = 0;
  for(int i = 0;i<nfiles;i++){
    string fstrp = base_fstrp;
    string fstrr = base_fstrr;
    std::ostringstream num;
    num<<i;
    if(nfiles > 1){
      fstrp = fstrp+"."+num.str();
      fstrr = fstrr+"."+num.str();
    }
    const char *fcharp = fstrp.c_str();
    const char *fcharr = fstrr.c_str();
    binifstream pfile(fstrp.c_str());
    binifstream rfile(fstrr.c_str());

    //confirm that files open
    if (pfile.fail()) {
      cerr<<"error: cannot open file '" <<fstrp<<"'"<<endl;
      exit (2030);
    }
    if (rfile.fail()) {
      cerr<<"error: cannot open file '" <<fstrr<<"'"<<endl;
      exit (2031);
    }


    //read the header information -- most of the information we don't care about
    cout<<"Reading particle data from file "<<fstrp<<endl;
    pfile>>buf; //it's f77_unformatted
    if (buf != 256){
      cout<<"Incorrect header format for file "<<fstrp<<", buf = "<<buf<<endl;
      return 0;
    }
    pfile.read((char *)(&header), sizeof(struct gadget_header));
    pfile>>buf;
    rfile>>buf>>buf>>buf>>buf;

    //read in the pos,vel,and rnn information
    float tmp_x,tmp_y,tmp_z, dist8;
#ifdef LONG64_IDS
    long int tid;
#else
    int tid;
#endif
    float vfac = sqrt(header.time);
    float xfac = 1./(sim.LengthUnit());

#ifdef DEBUG
    cout<<"Going to read "<<header.npart[1]<<" particles."<<endl;
#endif
    pfile>>buf;
    rfile>>buf;

    //define our structure to hold all particlese from this file
    struct triple{
      float x;
      float y;
      float z;
      float d8;
      float vx;
      float vy;
      float vz;
    } *tmp_pos;
    tmp_pos = new triple[header.npart[1]];
#ifdef LONG64_IDS
    long int *tmp_id;
    tmp_id = new long int[header.npart[1]];
#else
    int *tmp_id;
    tmp_id = new int[header.npart[1]];
#endif
    triple tmp_ent;

    //read in the particle positions and rdel values
    for (int ip=0;ip<header.npart[1];ip++){
      pfile>>tmp_x>>tmp_y>>tmp_z;
      tmp_pos[ip].x = tmp_x;
      tmp_pos[ip].y = tmp_y;
      tmp_pos[ip].z = tmp_z;
      rfile>>dist8;
      tmp_pos[ip].d8 = dist8;
    }
    pfile>>buf;

    //read in the particle positions
    pfile>>buf;
    for (int ip=0;ip<header.npart[1];ip++){
      pfile>>tmp_x>>tmp_y>>tmp_z;
      tmp_pos[ip].vx = tmp_x;
      tmp_pos[ip].vy = tmp_y;
      tmp_pos[ip].vz = tmp_z;
    }
    pfile>>buf;

    //read in the particle IDs -- need to be able to deal with long ids!!!!
    pfile>>buf;
    for(int ip=0;ip<header.npart[1];ip++){
      pfile>>tid;
      tmp_id[ip] = tid;
    }

    //read all the information -- closing our files
    pfile.close();
    rfile.close();

    //now loop through all read particles and save them into the correct classes
    for(int ip=0;ip<header.npart[1];ip++){ 
      tmp_x = tmp_pos[ip].x;
      tmp_y = tmp_pos[ip].y;
      tmp_z = tmp_pos[ip].z;
      Point xx(tmp_x*xfac,tmp_y*xfac,tmp_z*xfac);
#ifdef DEBUG
      if (ip == 0 || ip == 1)
	cout<<"First/Second particle position: "<<tmp_x<<" "<<tmp_y<<" "<<tmp_z<<endl;
#endif
      tmp_x = tmp_pos[ip].vx;
      tmp_y = tmp_pos[ip].vy;
      tmp_z = tmp_pos[ip].vz;
      Point vv(tmp_x*vfac,tmp_y*vfac,tmp_z*vfac);
      dist8 = tmp_pos[ip].d8;
      tid = tmp_id[ip];
      nr++;

      //Are we downsampling the particle distribution?
      if(drand48() > PARTICLE_FRACTION)
	continue;

      //We pretty much always do this -- catalot means run addgals!
#ifdef CATALOG
      Particle * particle = new Particle(xx,vv,dist8);
      //particle->Pid(nr);
      particle->Pid(tid);
//#ifdef DEBUG
      if(ip == 0 || ip == 1 || ip==2){
	cout<<"Info for first/second/third particle: "<<endl;
	xx.Print();
	vv.Print();
	cout<<"pid = "<<tid<<" "<<nr-1<<endl;
	cout<<"rnn = "<<dist8<<endl;
	cout<<particle->Ra()<<" "<<particle->Dec()<<" "<<particle->Zred()<<endl;
        cout<<"Distance from origin, redshift = "<<particle->R()<<" "<<particle->Zred()<<endl;
        cout<<"Ra, Dec = "<<particle->Ra()<<" "<<particle->Dec()<<endl;
      }
//#endif
      //Save or delete the particle
      if(particle->Save()){
	particles.push_back(particle);
	np++;
	if(particle->Zred()>zmax->GetVal())
	  zmax->SetVal(particle->Zred());
	if (np >= MAXSIZE){
	  cout<<"[ReadGadget]: ERROR!  Didn't allocate enough memory for particles!  MAXSIZE = "<<MAXSIZE<<endl;
	  exit(12);
	}

	if(particle->X() > maxpos[0])
	  maxpos[0] = particle->X();
	if(particle->Y() > maxpos[1])
	  maxpos[1] = particle->Y();
	if(particle->Z() > maxpos[2])
	  maxpos[2] = particle->Z();
	if(particle->X() < minpos[0])
	  minpos[0] = particle->X();
	if(particle->Y() < minpos[1])
	  minpos[1] = particle->Y();
	if(particle->Z() < minpos[2])
	  minpos[2] = particle->Z();

      }
      else{
	delete particle;
      }
#endif
#ifdef HALOASSIGN
      Particle * particle = new Particle(xx, vv, dist8);
      particles.push_back(particle);
      np++;
      if (np >= MAXSIZE)
	cout<<"[ReadGadget]: ERROR!  Didn't allocate enough memory for particles!  MAXSIZE = "<<MAXSIZE<<endl;
#endif
    }
    delete [] tmp_pos;
//#ifdef DEBUG
    cout<<" Total Particles read: "<<np<<endl;
    cout<<" Number of Particles in this file: "<<nr<<endl;
//#endif
  }
  cout<<"[read_gadget] "<<np<<" particles saved."<<endl;
  cout<<endl;

#ifdef DEBUG
  cout<<"min/max positions:"<<endl;
  cout<<"  x: "<<minpos[0]<<"/"<<maxpos[0]<<endl;
  cout<<"  x: "<<minpos[1]<<"/"<<maxpos[1]<<endl;
  cout<<"  x: "<<minpos[2]<<"/"<<maxpos[2]<<endl;
#endif

  return np;
}

//
// reads particles from file and returns vector of particle pointers
//
vector <Particle *> ReadGadgetParticles(int &nread){
  nread = 0;

  vector <Particle*> particles;
  // Reserve max number of particles to avoid memory difficulties
  particles.reserve(MAXSIZE);

  int nfiles = ReadGadgetnFiles();
  if (RAMAZ gt 180.0) then nfiles *= 2;

  cout<<"Reading "<<nfiles<<" gadget files."<<endl;

  nread = ReadGadget(nfiles,particles);

  cout<<"Done."<<endl;
  PRNT("[read_gadget] end:",nread);
  assert(particles.size()>0);

  return particles;

}
