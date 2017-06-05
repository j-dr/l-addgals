#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <list>
#include <strstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <functional>

#include "biniostream.h"
#include "point.h"
#include "particle.h"
#include "halo.h"
#include "singleton.h"
#include "functions.h"
#include "constants.h"
#include "gadget.h"
#include "ReadParameters.h"
#include "healpix_utils.h"

static const int MMAXINT=32767;

bool notInVolume(Particle *part)
{
  float tr = part->R();
  return ! ( ( RMIN_REAL <= tr ) && (tr <= RMAX_REAL) &&
	     ( RAMIN <= part->Ra() ) && ( part->Ra() < RAMAX ) &&
	     ( DECMIN <= part->Dec() ) && ( part->Dec() < DECMAX ) );
}

struct ToRing : public std::unary_function<long,long> {
  const int order_;
  ToRing(int o) : order_(o) {}
  long operator() (long pix) {return nest2ring(pix, order_);}
};

bool comparez(Particle *lhs, Particle *rhs)
{
  return lhs->Zred() < rhs->Zred();
}

bool comparera(Particle *lhs, Particle *rhs)
{
  return lhs->Ra() < rhs->Ra();
}

bool comparedec(Particle *lhs, Particle *rhs)
{
  return lhs->Dec() < rhs->Dec();
}

//
// reads all halos from file and returns vector of halo pointers
//

//routine for ROCKSTAR hlist files
vector <Halo*> ReadRockstarHalosHlist(void){
  vector <Halo*> halos;
  string filename, rnnfilename;

  cout<<"Reading Rockstar Halos (hlist format)."<<endl;
  filename = halofile;
  ifstream file(filename.c_str());
  if (file.fail()) {
    cerr<<"[ReadRockstarHalosHlist error]: cannot open "<<filename<<endl;    exit(1);
  }
  rnnfilename = rnn_halofile;
  ifstream rfile(rnnfilename.c_str());
  if (rfile.fail()) {
    cerr<<"[ReadRockstarHalosHlist error]: cannot open "<<rnnfilename<<endl;    exit(1);
  }

  string tmps;
  streampos oldpos;
  while(true){
  //  for(int i=0;i<43;i++){
    oldpos = file.tellg();  // stores the position
    getline(file, tmps);
    cout<<tmps<<endl;
    if (tmps[0]!='#') break;
  }
  file.seekg (oldpos);
  getline(rfile, tmps);
  cout<<tmps<<endl;

  while(true){
  //  for(int i=0;i<43;i++){
    oldpos = rfile.tellg();  // stores the position
    getline(rfile, tmps);
    cout<<tmps<<endl;
    if (tmps[0]!='0') break;
  }
  rfile.seekg (oldpos);

  int hcount = 0;
  int tcount = 0;

  while(file){
    int hid = 0;
    float xMpc,yMpc,zMpc,vx,vy,vz,mvir, rvir, vmax, vrms, rs, zred, rdel;
    float scale, dscale, sam_mvir, macc, mpeak, vacc, vpeak, last_mm, dummy;
    int ip, pid, did, nprog, upid, desc_pid, phantom, mmp, dummy_i;
    ip = 0;
    rdel=1;
    //file>>scale>>hid>>dscale>>did>>nprog>>pid>>upid>>desc_pid>>phantom>>m200c>>omvir>>r200c>>rs>>vrms>>mmp>>last_mm>>vmax>>xMpc>>yMpc>>zMpc>>vx>>vy>>vz>>macc>>mpeak>>vacc>>vpeak;
    getline(file, tmps);

    std::stringstream line_data(tmps);
    line_data>>scale>>hid>>dscale>>did>>nprog>>pid>>upid>>desc_pid>>phantom>>sam_mvir>>mvir>>rvir>>rs>>vrms>>mmp>>last_mm>>vmax>>xMpc>>yMpc>>zMpc>>vx>>vy>>vz;
    if (file.eof()) break;

    getline(rfile, tmps);
    sscanf(tmps.c_str(), "%d, %f", &hid, &rdel);
    //std::stringstream rline_data(tmps);
    //rline_data>>hid;
    //rline_data>>rdel;
    //cout<<hid<<rdel<<endl;
    //    if(pid >= 0) continue;
    //if(mvir < BCG_Mass_lim) continue;
    zred = cosmo.ZofR(sqrt(xMpc*xMpc + yMpc*yMpc + zMpc*zMpc));


    Point vel(vx,vy,vz);
    float xx, yy, zz;

    xx = xMpc;
    yy = yMpc;
    zz = zMpc;

    if (tcount<10)
      {
	cout << tmps << endl;
	cout << endl;
	cout << endl;
	cout << "xx, yy, zz, mvir, pid, hid, rvir, rdel: " << xMpc << " " << yMpc << " "<< zMpc << " "<< mvir << " "<< pid << " "<< hid << " "<< rvir << " "<< " " << rdel <<endl;
      }

    Point pos(xx, yy, zz);
    Halo * halo = new Halo(pos, vel, mvir, zred, ip, hid, vrms, rvir);
    halo->Dist8(rdel);
    if (pid < 0) pid = hid;
    halo->Host(pid);

#ifdef SNAPSHOT
    if(xx > 0 && yy > 0 && zz > 0 &&
       xx < sim.Boxsize()*BOXFR && yy < sim.Boxsize()*BOXFR &&
       zz < sim.Boxsize()*BOXFR)
      {
        halos.push_back(halo);
        hcount++;
      }
#else
    halos.push_back(halo);
    hcount++;
#endif
    tcount++;
  }
  file.close();
  cout<<"[ReadRockstarHalosHlist] Read "<<halos.size()<<" halos."<<endl;
  for (int i=0;i<10;i++) halos[i]->Print();
  return halos;
}


//routine for objects in ROCKSTAR ascii format
vector <Halo*> ReadRockstarHalos(void){
  vector <Halo*> halos;
  string filename, rnnfilename;

  cout<<"Reading Rockstar Halos."<<endl;
  cout<<"  halo file: "<<halofile<<endl;
  cout<<"  debugging."<<endl;
  filename = halofile;
  ifstream file(filename.c_str());
  if (file.fail()) {
    cerr<<"[ReadRockstarHalos error]: cannot open "<<filename<<endl;    exit(1);
  }
  rnnfilename = rnn_halofile;
  ifstream rfile(rnnfilename.c_str());
  if (rfile.fail()) {
    cerr<<"[ReadRockstarHalos error]: cannot open "<<rnnfilename<<endl;    exit(1);
  }

  //#  id, did, m200c, vmax, vrms, r200c, rs, np, x, y, z, vx, vy, vz, pid

  //read in first line of strings.
  string tmps;
  getline(file, tmps);
  cout<<tmps<<endl;
  getline(rfile, tmps);
  cout<<tmps<<endl;

  int hid = 0;
  int hcount = 0;
  while(file){
    double mvir;
    float xMpc,yMpc,zMpc,vx,vy,vz,m200c, r200c, vmax, vrms, rs, zred, rdel;
    int ip, hid, pid, did;
    file>>hid>>did>>m200c>>vmax>>vrms>>r200c>>rs>>ip>>xMpc>>yMpc>>zMpc>>vx>>vy>>vz>>pid;
    rfile>>hid>>rdel;
    zred = cosmo.ZofR(sqrt(xMpc*xMpc + yMpc*yMpc + zMpc*zMpc));

    //    if (hcount < 5){
    //      cout<<"Halo Position "<<hid<<" = "<<xMpc<<" "<<yMpc<<" "<<zMpc<<", rdel = "<<rdel<<endl;
    //    }

    if(m200c < BCG_Mass_lim) continue;
    if(pid >= 0) continue;
    Point vel(vx,vy,vz);
    float xx, yy, zz;

    xx = xMpc;
    yy = yMpc;
    zz = zMpc;

    Point pos(xx, yy, zz);
    //    if ((RMIN_REAL > pos.R()) || (pos.R() > RMAX_REAL)) continue;

    Halo * halo = new Halo(pos, vel, m200c, zred, ip, hid, vrms, r200c);
    halo->Host(pid);
    halo->Dist8(rdel);
    //halo->Particle(0);

#ifdef SNAPSHOT
    if(xx > 0 && yy > 0 && zz > 0 &&
       xx < sim.Boxsize()*BOXFR && yy < sim.Boxsize()*BOXFR &&
       zz < sim.Boxsize()*BOXFR)
      {
        halos.push_back(halo);
        hcount++;
      }
#else
    halos.push_back(halo);
    hcount++;
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

      //if(mvir < BCG_Mass_lim) continue;
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

  sim.ParticleMass(header.mass[1]*1e10);

  pfile.close();
  sim_redshift = header.redshift;
#ifdef SNAPSHOT
  ZREDMIN = sim_redshift;
  ZREDMAX = sim_redshift;
#endif
  return header.num_files;
}

int ReadGadget(int nfiles, vector <Particle *> &particles){
  gadget_header header;
  int buf;
  long int nr = 0;
  int np = 0;
  double t1, t2, t3, TimeRead;

  //define out input files -- one for the particles, one for rnn
  string base_fstrp = simulationfile;
  string base_fstrr = rnnfile;
  string base_fstrh;
  if (read_hod) base_fstrh = hodfile;

  float maxpos[3];
  float minpos[3];
  maxpos[0] = maxpos[1] = maxpos[2] = 0.;
  minpos[0] = minpos[1] = minpos[2] = 10000000.;

  nr = 0;
  for(int i = 0;i<nfiles;i++){
    t1 = clock();
    string fstrp = base_fstrp;
    string fstrr = base_fstrr;
    string fstrh = base_fstrh;
    std::ostringstream num;
    num<<i;
    if(nfiles > 1){
      fstrp = fstrp+"."+num.str();
      fstrr = fstrr+"."+num.str();
      fstrh = fstrh+"."+num.str();
    }
    const char *fcharp = fstrp.c_str();
    const char *fcharr = fstrr.c_str();
    const char *fcharh = fstrh.c_str();
    std::ifstream pfile(fstrp.c_str());
    std::ifstream rfile(fstrr.c_str());
    ifstream hfile(fstrh.c_str());

    //confirm that files open
    if (pfile.fail()) {
      cerr<<"error: cannot open file '" <<fstrp<<"'"<<endl;
      exit (2030);
    }
    if (rfile.fail()) {
      cerr<<"error: cannot open file '" <<fstrr<<"'"<<endl;
      exit (2031);
    }
    if (read_hod) {
      if (hfile.fail()) {
        cerr<<"error: cannot open file '" <<fstrh<<"'"<<endl;
        exit (2032);
      }
    }

    //read the header information -- most of the information we don't care about
    cout<<"Reading particle data from file "<<fstrp<<endl;
    pfile.read((char*) &buf, 4); //it's f77_unformatted
    if (buf != 256){
      cout<<"Incorrect header format for file "<<fstrp<<", buf = "<<buf<<endl;
      return 0;
    }

#ifdef DEBUG
    cout<<"Reading header of file"<<endl;
#endif

    pfile.read((char *)(&header), sizeof(struct gadget_header));
    pfile.read((char*) &buf, 4);
    rfile.read((char*) &buf, 4);
    rfile.read((char*) &buf, 4);
    rfile.read((char*) &buf, 4);
    rfile.read((char*) &buf, 4);
    string tmps;
    if (read_hod) {
      //hfile>>tmps;
      getline(hfile,tmps);
      //cout<<tmps<<endl;
    }

    //read in the pos,vel,and rnn information
#ifdef DEBUG
    cout<<"Reading in the pos, vel and rnn info"<<endl;
#endif
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
    pfile.read((char*) &buf, 4);
    rfile.read((char*) &buf, 4);
#ifdef DEBUG
    cout<<"creating structures"<<endl;
#endif
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
    struct xyz{
      float x;
      float y;
      float z;
    } *tmp_cart;
    tmp_cart = new xyz[header.npart[1]];
    float *tmp_d8;
    tmp_d8 = new float[header.npart[1]];
#ifdef LONG64_IDS
    long int *tmp_id;
    tmp_id = new long int[header.npart[1]];
    int ID_SIZE = sizeof(long int);
#else
    int *tmp_id;
    tmp_id = new int[header.npart[1]];
    int ID_SIZE = sizeof(int);
#endif
    triple tmp_ent;

    //define some structures for storing the particle HOD information
    long int *tmp_phid;
    float *tmp_pmvir;
    float *tmp_prhalo;
    tmp_phid = new long int[header.npart[1]];
    tmp_pmvir = new float[header.npart[1]];
    tmp_prhalo = new float[header.npart[1]];
#ifdef DEBUG
    cout<<"read in the particle positions and rdel values"<<endl;
#endif
    //read in the particle positions and rdel values
    pfile.read((char*) &tmp_cart[0], sizeof(float)*3*header.npart[1]);
    rfile.read((char*) &tmp_d8[0], sizeof(float)*header.npart[1]);
    cout<<"particle reading done"<<endl;
    for (int ip=0;ip<header.npart[1];ip++){
      tmp_pos[ip].x = tmp_cart[ip].x;
      tmp_pos[ip].y = tmp_cart[ip].y;
      tmp_pos[ip].z = tmp_cart[ip].z;
      tmp_pos[ip].d8 = tmp_d8[ip];
      if (tmp_d8[ip] < 0 || tmp_d8[ip] > 22){
	cout<<"Error with read in rnn = "<<tmp_d8[ip]<<" for particle "<<ip<<endl;
      }
      if (read_hod){
	int thid;
	float tmvir, trhalo;
	hfile>>thid>>tmvir>>trhalo;
	tmp_phid[ip] = thid;
	tmp_pmvir[ip] = tmvir;
	tmp_prhalo[ip] = trhalo;
	//if (ip < 10) cout<<thid<<" "<<tmvir<<" "<<trhalo<<endl;
      }
    }
#ifdef DEBUG
    cout<<"Particle positions being read"<<endl;
#endif
    pfile.read((char*) &buf, 4);

    //read in the particle positions
    pfile.read((char*) &buf, 4);
    pfile.read((char*) &tmp_cart[0], sizeof(float)*3*header.npart[1]);
    for (int ip=0;ip<header.npart[1];ip++){
      tmp_pos[ip].vx = tmp_cart[ip].x;
      tmp_pos[ip].vy = tmp_cart[ip].y;
      tmp_pos[ip].vz = tmp_cart[ip].z;
    }
    pfile.read((char*) &buf, 4);

    //read in the particle IDs -- need to be able to deal with long ids!!!!
    pfile.read((char*) &buf, 4);
    pfile.read((char*) &tmp_id[0], header.npart[1]*ID_SIZE);

#ifdef DEBUG
    cout<<"closing files"<<endl;
#endif
    //read all the information -- closing our files
    pfile.close();
    rfile.close();
    if (read_hod) hfile.close();

    t2 = clock();

/*
    //initialize a dummy particle
    tmp_x = tmp_pos[0].x;
    tmp_y = tmp_pos[0].y;
    tmp_z = tmp_pos[0].z;
    Point xx(tmp_x*xfac,tmp_y*xfac,tmp_z*xfac);
    tmp_x = tmp_pos[0].vx;
    tmp_y = tmp_pos[0].vy;
    tmp_z = tmp_pos[0].vz;
    Point vv(tmp_x*vfac,tmp_y*vfac,tmp_z*vfac);
    dist8 = tmp_pos[0].d8;
    Particle * particle = new Particle(xx,vv,dist8);
*/
    //now loop through all read particles and save them into the correct classes
    for(int ip=0;ip<header.npart[1];ip++){
      tmp_x = tmp_pos[ip].x;
      tmp_y = tmp_pos[ip].y;
      tmp_z = tmp_pos[ip].z;
      Point xx(tmp_x*xfac,tmp_y*xfac,tmp_z*xfac);
//#ifdef DEBUG
      if (ip == 0 || ip == 1 || ip == header.npart[i]-1)
	cout<<"First/Second/Last particle position, rnn: "<<tmp_x<<" "<<tmp_y<<" "<<tmp_z<<" "<<tmp_pos[ip].d8<<endl;
//#endif
      tmp_x = tmp_pos[ip].vx;
      tmp_y = tmp_pos[ip].vy;
      tmp_z = tmp_pos[ip].vz;
      Point vv(tmp_x*vfac,tmp_y*vfac,tmp_z*vfac);
      dist8 = tmp_pos[ip].d8;
      tid = tmp_id[ip];
/*
      if (!(dist8 < 15.0))
      {
	cout<<"Error in dist8!  Particle "<<ip<<" in file "<<i<<" !!!"<<endl;
	cout<<"Position: "<<tmp_pos[ip].x<<" "<<tmp_pos[ip].y<<" "<<tmp_pos[ip].z<<endl;
	cout<<"Density: "<<dist8<<endl;
	cout<<"  (previous particle had density "<<tmp_pos[ip-1].d8<<")"<<endl;
      }
*/
      nr++;

      //Are we downsampling the particle distribution?
      if(drand48() > PARTICLE_FRACTION)
	continue;

      //we are are now checking this in assignment
      //is the particle above the mhost cut?
      //if(read_hod)
	//if (tmp_pmvir[ip] > mhost_cut && dist8 < rnn_cut)
	  //continue;

      //We pretty much always do this -- catalog means run addgals!
#ifdef CATALOG
      Particle * particle = new Particle(xx,vv,dist8);
      particle->Pid(tid);
      if (read_hod) {
	particle->Hid(tmp_phid[ip]);
	particle->Mhost(tmp_pmvir[ip]);
	//cout<<tmp_phid[ip]<<" "<<tmp_pmvir[ip]<<endl;
      }
/*
      particle->PosAssign(xx);
      particle->VelAssign(vv);
      particle->DensAssign(dist8);
      particle->Pid(tid);
*/
#ifdef DEBUG
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
#endif
      //Save or delete the particle
      if(particle->Save()){
        //if (particle->Zred() > ZREDMAX) cout<<"Saved particle redshift out of range: "<<particle->Zred()<<endl;
	particles.push_back(particle);
	np++;
#ifdef SNAPSHOT
	if(particle->Z()>zmax->GetVal())
	  zmax->SetVal(particle->Z());
#else
	if(particle->Zred()>zmax->GetVal())
	  zmax->SetVal(particle->Zred());
#endif
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
    delete [] tmp_cart;
    delete [] tmp_d8;
    delete [] tmp_id;
    t3 = clock();
//#ifdef DEBUG
    cout<<" Total Particles saved: "<<np<<endl;
    cout<<" Number of Particles in this file: "<<nr<<endl;
    TimeRead = (t2-t1)/CLOCKS_PER_SEC;
    cout<<" Time to read the file: "<<TimeRead<<endl;
    TimeRead = (t3-t2)/CLOCKS_PER_SEC;
    cout<<" Time to save the particle: "<<TimeRead<<endl;
    TimeRead = (t3-t1)/CLOCKS_PER_SEC;
    cout<<" Time to process the file: "<<TimeRead<<endl;
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

#ifdef DEBUGLC
  string pfname = "lcdata.dat";
  ofstream pfile(pfname.c_str());
  for (vector<Particle *>::iterator pitr=particles.begin(); pitr!=particles.end(); pitr++)
    {
      Particle *p = *pitr;
      p->Write(pfile);//pfile << pitr->X() << " " << pitr->Y() << " " pitr->Z() << " " pitr->VX() << " " << pitr->VY() << " " pitr->VZ()
      //   << " " << pitr->Zred() << endl;
    }
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
#ifdef BCC
  if (RAMAX > 90.0 && nfiles > 2) nfiles *= 2;
#endif

  cout<<"Reading "<<nfiles<<" gadget files."<<endl;

  nread = ReadGadget(nfiles,particles);

  cout<<"Done."<<endl;
  PRNT("[read_gadget] end:",nread);
  sim.Np(nread);
  assert(particles.size()>0);

  return particles;

}

#ifdef BCC
struct io_header readLCCellHeader(std::string fname)
{
  std::ifstream pfile(fname.c_str());
  struct io_header header;

  if (pfile.fail()) {
    cerr<<"error: cannot open file '" <<fname<<"'"<<endl;
    exit(2031);
  }

  pfile.read((char *)(&header), sizeof(struct io_header));

  return header;
}

unsigned int getIndexNSide()
{
  int r = 0;
  struct io_header header;

  while (true) {
    std::ostringstream convert;
    convert << datadir << simlabel << "_" << r << "_0";
    std::string fname = convert.str();
    std::ifstream pfile(fname.c_str());

    if (pfile.fail()) {
      r++;
      continue;
    }

    pfile.read((char *)(&header), sizeof(struct io_header));
    return header.nside;
  }
}

unsigned int getFileNSide()
{
  int r = 0;
  struct io_header header;

  while (true) {
    std::ostringstream convert;
    convert << datadir << simlabel << "_" << r << "_0";
    std::string fname = convert.str();
    std::ifstream pfile(fname.c_str());

    if (pfile.fail()) {
      r++;
      continue;
    }

  pfile.read((char *)(&header), sizeof(struct io_header));
  return header.filenside;
  }
}

void getRadialBins(int &minrb, int &maxrb)
{
  float r_zmin = cosmo.RofZ(ZREDMIN);
  float r_zmax = cosmo.RofZ(ZREDMAX);
  minrb = r_zmin/25;
  maxrb = r_zmax/25;
}

std::vector<long> getFilePixels(string datadir, string simlabel, long pix,
				int r, long order1_)
{
  /*
    inputs:
    pix -- The pixel that we are adding galaxies to
    order1_ -- The healpix order that addgals is using
    order2_ -- The healpix order that the LC radial bin uses
  */

  string fname;
  struct io_header header;
  int temp;
  long order2_=0;
  long pc = 0;

  //Look for first pixel with particles in it for this radial bin
  //get the file nside from this file
  std::ostringstream convert;
  convert << datadir << simlabel << "_" << r << "_" <<pc;
  fname = convert.str();
  header = readLCCellHeader(fname);

  temp = header.filenside;
  while(temp >>= 1) order2_++;

  pix = ring2nest(pix, order1_);

#ifdef DEBUG_PIXLC
  cout << "order1_, order2_: " << order1_ << order2_ << endl;
  cout << "Working cell nested pixel number: " << pix << endl;
#endif


  std::vector<long> fpix;
  if (order2_ < order1_)
    {
      //Only one file contains our pixel
      fpix.push_back(lower_nest(pix, order1_, order2_));
    }
  else if (order2_ > order1_)
    {
      fpix.resize(1 << ( 2 * ( order2_ - order1_ ) ) );
      higher_nest(pix, order1_, order2_, &fpix[0]);
    }
  else
    {
      fpix.push_back(pix);
    }

  return fpix;
}

long getNparts(int minrb, int maxrb, long order1_, long order2_, vector<long> &pidx)
{
  int i,r,count,temp;
  long np,j;
  long pc;
  int npix = 12*(1<<(2*order2_));
  long fnpix=9999999999;
  long forder;
  long nparts = 0;
  vector<long> idx(npix);
  std::string fname;
  std::vector<long> fpix;
  struct io_header header;

  // Only using first two octants
  for (r=minrb; r<=maxrb; r++)
    {
      //Determine which file pixels contain particles in the
      //healpix cell we are using now
      fpix = getFilePixels(datadir, simlabel, PixelNum, r, order1_);

      //Add up parts from each lightcone pixel associated with PixelNum
      for (vector<long>::iterator itr=fpix.begin(); itr!=fpix.end(); itr++)
	{
	  pc = *itr;
	  std::ostringstream convert;
	  convert << datadir << simlabel << "_" << r << "_" << pc;
	  fname = convert.str();
	  std::ifstream pfile(fname.c_str());

	  //Sometimes a pixel will be empty, if so raise a warning but continue
	  if (pfile.fail()) {
	    cerr<<"warning: cannot open file '" <<fname<<"'"<<endl;
	    continue;
	  }

	  pfile.read((char *)(&header), sizeof(struct io_header));
	  if (header.npart == 0) continue;

	  //Read in the index, make sure it agrees with npart from header
	  pfile.read((char *)(&idx[0]), sizeof(long)*12*(1<<(2*order2_)));
	  partial_sum(idx.begin(), idx.end(), idx.begin());
	  assert(header.npart == idx[idx.size()-1]);

	  long fpos = 0;
	  //Add up the particles that are in the right peano indices
	  for (vector<long>::iterator itr=pidx.begin(); itr!=pidx.end(); itr++)
	    {
	      if (*itr==0)
		{
		  np = idx[*itr];
		}
	      else
		{
		  np = idx[*itr] - idx[(*itr)-1];
		}
	      nparts += np;
	    }
	}
    }
  return nparts;
}


vector <Particle *> ReadGadgetLCCell()
{
  int minrb,maxrb;
  int i,j,r,buf;
  int indexnside, temp;
  int nbins;
  long order1_ = 0, order2_ = 0;
  long tid, accum;
  long step, pnp, nparts;
  float td, xfac, vfac;
  float tstart, tend;
  long pix=0;
  long fnpix=9999999999;
  long fnp=0;
  std::vector<long> fpix;
  std::string fname;
  std::string rnnfname, hinfofname;
  struct io_header header;

  struct triple{
    float x;
    float y;
    float z;
  } cart;

  hrow hdata;

  cout << "Determining radial bins" << endl;
  getRadialBins(minrb, maxrb);
  cout << "Radial bins " << minrb << " " << maxrb << endl;
  indexnside = getIndexNSide();
  cout << "IndexNside: " << indexnside << endl;
  nbins = maxrb-minrb;

  temp = nSide;
  while(temp >>= 1) order1_++;

  temp = indexnside;
  while(temp >>= 1) order2_++;

#ifdef DEBUG_PIXLC
  cout << "order1_, order2_: " << order1_ << " " << order2_ << endl;
#endif

  vector<long> idx(12*(1<<(2*order2_)));
  vector<long> pidx(1<<(2*(order2_-order1_)));

  ring2peanoindex(PixelNum, order1_, order2_, pidx);
  nparts = getNparts(minrb, maxrb, order1_, order2_, pidx);
  cout << "[ReadGadgetLCCell] Reading in " << nparts << "particles" <<endl;
  vector<Particle *> parts(nparts);
  tstart = clock();
  // Only using first two octants

  //vfac = sqrt(header.time); does this do anything with LCs?
  xfac = 1./(sim.LengthUnit());
  accum = 0;
  for (r=minrb; r<=maxrb; r++)
    {


      fpix = getFilePixels(datadir, simlabel, PixelNum, r, order1_);


      for (vector<long>::iterator itr=fpix.begin(); itr!=fpix.end(); itr++)
	{
	  pix = *itr;
#ifdef DEBUG_PIXLC
	  cout << "Reading file: " << datadir << simlabel << "_" << r
	       << "_" << pix << endl;
#endif
	  std::ostringstream convert;
	  convert << datadir << simlabel << "_" << r
		  << "_" << pix;
	  fname = convert.str();
	  std::ifstream pfile(fname.c_str());
	  if (pfile.fail()) {
	    cerr<<"error: cannot open file '" <<fname<<"'"<<endl;
	    exit(3031);
	  }

	  pfile.read((char *)(&header), sizeof(struct io_header));
	  if (header.npart == 0) continue;

	  std::ostringstream rconvert;
	  rconvert << datadir << "rnn_" << simlabel << "_" << r
		   << "_" << pix;
	  rnnfname = rconvert.str();
	  std::ifstream rfile(rnnfname.c_str());
	  if (rfile.fail()) {
	    cerr<<"error: cannot open file '" <<rnnfname<<"'"<<endl;
	    exit (2031);
	  }

      std::ostringstream hconvert;
	  hconvert << datadir << "hinfo_" << simlabel << "_" << r
		   << "_" << pix;
	  hinfofname = hconvert.str();
	  std::ifstream hfile(hinfofname.c_str());
	  if (hfile.fail()) {
	    cerr<<"error: cannot open file '" <<hinfofname<<"'"<<endl;
	    exit (2031);
	  }


	  pfile.read((char *)(&idx[0]), sizeof(long)*12*(1<<(2*order2_)));
	  partial_sum(idx.begin(), idx.end(), idx.begin());
	  assert(header.npart == idx[idx.size()-1]);
	  cout << "number of particles in this file is " << header.npart << endl;


	  rfile.read((char*) &buf, 4);
	  rfile.read((char*) &buf, 4);
	  rfile.read((char*) &buf, 4);
	  rfile.read((char*) &buf, 4);
	  rfile.read((char*) &buf, 4);
	  cout << "number of densities in this file is " << buf/sizeof(float) << endl;

	  // read in positions and densities
	  long fpos = 0;
	  fnp = 0;
	  for (vector<long>::iterator itr=pidx.begin(); itr!=pidx.end(); itr++)
	    {
	      if (*itr==0)
		{
		  step = 0;
		  pnp = idx[*itr];
		}
	      else if (itr==pidx.begin())
		{
		  step = idx[(*itr)-1];
		  pnp = idx[*itr] - idx[(*itr)-1];
		}
	      else
		{
		  step = idx[(*itr)-1] - idx[*(itr-1)];
		  pnp = idx[*itr] - idx[(*itr)-1];
		}

	      fpos += step + pnp;
	      pfile.seekg( 3*step*sizeof(float), pfile.cur );
	      rfile.seekg( step*sizeof(float), rfile.cur );
	      hfile.seekg( step*(sizeof(hrow)), hfile.cur);
	      for (j=0; j<pnp; j++)
		{
		  //cout << j << endl;
		  pfile.read((char *)&cart, sizeof(struct triple));
		  rfile.read((char *)&td, sizeof(float));
		  hfile.read((char *)&hdata, sizeof(struct hrow));
		  Point xx(cart.x*xfac,cart.y*xfac,cart.z*xfac);
		  parts[accum+fnp] = new Particle();
		  parts[accum+fnp]->PosAssign(xx);
		  parts[accum+fnp]->DensAssign(td);
		  parts[accum+fnp]->HaloAssign(hdata);
		  fnp+=1;
		}
	    }
	  //Seek to the end of the positions
	  step = (idx.back())-idx[pidx.back()];
	  pfile.seekg( 3*step*sizeof(float), pfile.cur );
	  cout << "Read in " << fnp << " particles" << endl;
	  // read in velocities
	  fnp = 0;
	  for (vector<long>::iterator itr=pidx.begin(); itr!=pidx.end(); itr++)
	    {
	      if (*itr==0)
		{
		  step = 0;
		  pnp = idx[*itr];
		}
	      else if (itr==pidx.begin())
		{
		  step = idx[(*itr)-1];
		  pnp = idx[*itr] - idx[(*itr)-1];
		}
	      else
		{
		  step = idx[(*itr)-1] - idx[*(itr-1)];
		  pnp = idx[*itr] - idx[(*itr)-1];
		}

	      pfile.seekg( 3*step*sizeof(float), pfile.cur );
	      for (j=0; j<pnp; j++)
		{
		  pfile.read((char *)&cart, sizeof(struct triple));
		  Point vv(cart.x,cart.y,cart.z);
		  parts[accum+fnp]->VelAssign(vv);
		  fnp+=1;
		}
	    }
	  //Seek to the end of the velocities
	  step = idx.back()-idx[pidx.back()];
	  pfile.seekg( 3*step*sizeof(float), pfile.cur );
	  // read in particle ids
	  fnp = 0;
	  for (vector<long>::iterator itr=pidx.begin(); itr!=pidx.end(); itr++)
	    {
	      if (*itr==0)
		{
		  step = 0;
		  pnp = idx[*itr];
		}
	      else if (itr==pidx.begin())
		{
		  step = idx[(*itr)-1];
		  pnp = idx[*itr] - idx[(*itr)-1];
		}
	      else
		{
		  step = idx[(*itr)-1] - idx[*(itr-1)];
		  pnp = idx[*itr] - idx[(*itr)-1];
		}

	      pfile.seekg( step*sizeof(long), pfile.cur );
	      for (j=0; j<pnp; j++)
		{
		  pfile.read((char *)&tid, sizeof(long));
		  parts[accum+fnp]->Pid(tid);
		  fnp+=1;
		}
	    }
	  accum += fnp;
	  cout << "Have read " << accum << " particles in total" << endl;
	  pfile.close();
	  rfile.close();
	}
    }
  assert( accum == nparts );
#ifdef DEBUG_PIXLC
  vector<Particle*>::iterator maxzpp = std::max_element(parts.begin(), parts.end(), comparez);
  vector<Particle*>::iterator maxrap = std::max_element(parts.begin(), parts.end(), comparera);
  vector<Particle*>::iterator maxdecp = std::max_element(parts.begin(), parts.end(), comparedec);
  cout << "Maximum particle redshift, ra, dec before cut is : " << parts[maxzpp-parts.begin()]->Zred()
       << " " << parts[maxrap-parts.begin()]->Ra() << " " << parts[maxdecp-parts.begin()]->Dec() << endl;
  vector<Particle*>::iterator minzpp = std::min_element(parts.begin(), parts.end(), comparez);
  vector<Particle*>::iterator minrap = std::min_element(parts.begin(), parts.end(), comparera);
  vector<Particle*>::iterator mindecp = std::min_element(parts.begin(), parts.end(), comparedec);
  cout << "Minimum particle redshift, ra, dec before cut is : " << parts[minzpp-parts.begin()]->Zred()
       << " " << parts[minrap-parts.begin()]->Ra() << " " << parts[mindecp-parts.begin()]->Dec() << endl;
  cout << "Cut quantities are RAMIN, RAMAX, DECMIN, DECMAX: " << RAMIN << " " << RAMAX << " " << DECMIN
       << " " << DECMAX << endl;
  
#endif
  //Get rid of particles that fall outside of redshift range
  parts.erase( std::remove_if( parts.begin(), parts.end(), notInVolume ), parts.end());

  //Get maximum redshift
  vector<Particle*>::iterator maxzp = std::max_element(parts.begin(), parts.end(), comparez);

  cout << "Maximum particle redshift is : " << parts[maxzp-parts.begin()]->Zred() << endl;
  zmax->SetVal(parts[maxzp-parts.begin()]->Zred());

  tend = clock();

  cout << "[ReadGadgetLCCell] Reading particles took " << (tend-tstart)/CLOCKS_PER_SEC << " seconds." <<endl;
  
  cout << "Halo information from first few particles: " << endl;
  for (i=0; i<5; i++)
    {
      cout << "Hid, RHalo, RVir, MVir : " << parts[i]->Hid() << " " << parts[i]->RHalo() << " " << parts[i]->RVir() << " " << parts[i]->MVir() << endl;
    }

#ifdef DEBUGLC
  string pfname = "lcdata.dat";
  ofstream pfile(pfname.c_str());
  for (vector<Particle *>::iterator pitr=parts.begin(); pitr!=parts.end(); pitr++)
    {
      Particle *p = *pitr;
      p->Write(pfile);//pfile << pitr->X() << " " << pitr->Y() << " " pitr->Z() << " " pitr->VX() << " " << pitr->VY() << " " pitr->VZ()
      //   << " " << pitr->Zred() << endl;
    }
#endif

  return parts;
}
#endif
