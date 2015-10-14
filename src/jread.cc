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

static const int MMAXINT=32767;

ofstream errfile("errors.dat");


//
// Converts short int position into box length units
//
float PositionConvert(short int x, int cube, int xs){
  float pos;
  if(simulation==HV){
    pos = (((float) x)/((float) MMAXINT)+((float) cube))/((float) BOXES);
    
    if(xs<=7){
	//out<i<<" "<<j<<" "<<k<<" "<<pos<<" "<<;
	//cout<<cube<<" "<<xs<<" "<<pos<<" ";
      pos = 0.5-pos;
      //cout<<pos<<endl;
      }
      else{
	//    cout<<"2 "<<(box==HVPOW)<<" "<<cube<<" "<<xs<<" "<<pos<<" ";
	pos = pos-0.5;
	//cout<<pos<<endl;
      }
    }
  return pos;
}

//
// reads all halos from file and returns vector of halo pointers
//
vector <Halo*> ReadHVHalos(void){
  vector <Halo*> halos;

  string filename;
#ifdef HVL
  if (box == HVLIGHT)
    filename = Dir()+"lcdm.MS.msort12";
  else if (box == HVZ0)
    filename = Dir()+"lcdm.z0.msort12";
  else if (box == HVPOW)
    filename = Dir()+"lcdm.POW.msort22";
  else{
    cerr<<"[read_cube error]: HV box type not defined"<<box<<endl;    exit(1);
  }
#endif
  ifstream file(filename.c_str());
  if (file.fail()) {
    cerr<<"[read_cube error]: cannot open "<<filename<<endl;    exit(1);
  }
  //#  m15    zredsigma ip    xGpc    yGpc    zGpc    vx    vy    vzsiglos rdelta
  //#

  //read in first two lines of strings.
  for(int i=0; i<13;i++){
  string tmps;
  file>>tmps;
  }
  int hid = 0;
  while(file){
    double m15;
    float zred,sigma,xGpc,yGpc,zGpc,vx,vy,vz,siglos,rdelta;
    int ip, pid; 
    file>>m15>>zred>>sigma>>ip>>xGpc>>yGpc>>zGpc>>vx>>vy>>vz>>siglos>>rdelta>>pid;
    
    Point vel(vx,vy,vz);   
    float xx, yy, zz;
    
    //cout<<"you better put this back!<<endl;
    //    if(normalization == SHIFT){
    // float offset_cells = xstart-8;
    // xx = xGpc*1000.-1500.0-offset_cells*187.5;
    // yy = yGpc*1000.-1500.0-offset_cells*187.5;
    // zz = zGpc*1000.-1500.0-offset_cells*187.5;
      //	cout<<xGpc*1000-1500.0<<" "<<xx<<endl;
    //}
    //else{
    if(box==HVPOW){
      xx = xGpc*1000.f;
      yy = yGpc*1000.f;
      zz = zGpc*1000.f;
    }
    else{
    xx = xGpc*1000.f-1500.0f;
    yy = yGpc*1000.f-1500.0f;
    zz = zGpc*1000.f-1500.0f;
    }
    //}
    
    Point pos(xx, yy, zz);
    Halo * halo = new Halo(pos, vel, m15*1e15, zred, ip, hid, siglos, rdelta);
    
    //  if(ip > 0){  this was a bug introducer...
    // cout<<m15<<" "<<pos.X()<<" "<<vx<<" "<<hid<<endl;

    halos.push_back(halo);
    hid++;
      
      //}
  }
  file.close();
  cout<<"[read_cube] Read "<<halos.size()<<" halos."<<endl;
  return halos;
}


//
// reads particles from one data cube (i,j,k)
//
int ReadCube(int i, int j, int k, vector <Particle *> &particles, 
	      int xstart, int ystart, int zstart){
  //cout<<"[read_cube] looking for cube "<<i<<" "<<j<<" "<<k<<endl; 
  //string dir = Dir();
  string pdir = Dir()+"cubes/position/";
  string vdir = Dir()+"cubes/velocity/";
  string rdir = Dir()+"rnn/";
  string hdir = Dir()+"rnn/";
  
 string fstrp = pdir+"p."+MakeString(i,2)+"."
   +MakeString(j,2)+"."+MakeString(k,2);
 string fstrv = vdir+"v."+MakeString(i,2)+"."
   +MakeString(j,2)+"."+MakeString(k,2);
 string fstrh = rdir+"h."+MakeString(i,2)+"."
   +MakeString(j,2)+"."+MakeString(k,2);
 string fstrhid = hdir+"ahid."+MakeString(i,2)+"."
   +MakeString(j,2)+"."+MakeString(k,2);

 //cout <<dir<<" "<<fstrp<<" "
 //    <<fstrv<<" "
 //    <<fstrh<<" "
 //    <<fstrhid<<" "<<endl;

  binifstream pfile(fstrp.c_str());
  binifstream vfile(fstrv.c_str());
  binifstream hfile(fstrh.c_str());
  ifstream hidfile(fstrhid.c_str());


  //confirm that files open
  if (pfile.fail()) {
    cerr<<"error: cannot open file '" <<fstrp<<"'"<<endl;
    return 0;
  }
  //else cerr<<"reading file '"<<fstrp<<"'"<<endl;
  if (vfile.fail()) {
    cerr<<"error: cannot open file '" <<fstrv<<"'"<<endl;
    return 0;
  }
  //else cerr<<"reading file '"<<fstrv<<"'"<<endl;
  
  if (hfile.fail()) {
    cerr<<"error: cannot open file '" <<fstrh<<"'"<<endl;
    return 0;
  }
  //else cout<<"reading file '"<<fstrh<<"'"<<endl;

#ifdef DEBUG
  cout<<"[read_cube] reading files "<<fstrp<<" "<<fstrv<<" "<<fstrh<<" "<<fstrhid<<endl;
#endif

  //#ifndef SPARC
  //convert endian-ness on intel & dec machines
  //you need this on ghostwheel.
  //ghostwheel is little endian.  the files are big endian.
  //therefore you need to read them as such.
  
#ifdef LITTLEENDIAN
  if(box!=HVPOW){
    pfile.bigendian();
    vfile.bigendian();
    hfile.bigendian();
  }
#endif
  
  //#endif

  short int tmp_x, tmp_y, tmp_z, tmp_dens;
  unsigned int np = 0;
  unsigned int nr = 0;

  int hid;

  while((pfile)&&(vfile)&&(hfile)){
  //while((pfile)&&(vfile)){
    //** Read 3 short ints from particle file
    //** Convert to box positions as they are passed to constructor
    pfile>>tmp_x>>tmp_y>>tmp_z;
    Point xx(PositionConvert(tmp_x,i,xstart), 
	     PositionConvert(tmp_y,j,ystart), 
	     PositionConvert(tmp_z,k,zstart));

    //** Read 3 short ints from velocity file, save as ShortPoint
    vfile>>tmp_x>>tmp_y>>tmp_z;
    Point vv(tmp_x, tmp_y, tmp_z);

    //** Read 1 short int from density file
    hfile>>tmp_dens;
    nr++;

    float hminabs = 1. / (16.*2.*MMAXINT);
    float sfaci =  (log(1.) - log(hminabs))/  (2.*MMAXINT);
    float dist = hminabs * exp( sfaci* (tmp_dens +MMAXINT) );
    float dist8= dist*sim.LengthUnit();
    hidfile>>hid;   


    Particle * particle = new Particle(xx, vv, dist8);
    if(hid>=0){
      particle->Hid(hid);
    }
    if(particle->Save()){
      particles.push_back(particle);
      np++;
      if(particle->Zred()>zmax->GetVal())
	  zmax->SetVal(particle->Zred());
    }
    else{
      delete particle;
    }
  }
  if((pfile)||(vfile)||(hfile)){
    cerr<<"[read_cube] Problem: files are different lengths"<<pfile<<" "<<vfile<<" "<<hfile<<endl;
    int nextra = 0;
    while(hfile){
      float tt;
      hfile>>tt;
      nextra++;
    }
    PRNT("read_cube",nextra);
    // cout<<"ERROR"
    errfile<<i<<" "<<j<<" "<<k<<" "<<nr<<" "<<nextra<<endl;
    cout<<"[read_cube] "<<nr<<" particles read in this cube"<<"...;"
	<<pfile<<" "<<vfile<<" "<<hfile<<endl;
    return 0;
    
  }
  cout<<"[read_cube] "<<i<<" "<<j<<" "<<k<<" "<<nr<<" particles read in this cube"<<"...";
  cout<<"[read_cube] "<<np<<" particles saved."<<endl;
  //  cout<<particles.size()<<" total particles read"<<endl;

  //cerr<<"  *****  boundaries:"<<i<<" "<<j<<" "<<i/16.<<" "<<(i+1)/16.<<" "
  //  <<minx<<" "<<maxx<<endl;
  PRNTV(nr);
 // exit(0);
  return nr;
}

//
// reads particles from file and returns vector of particle pointers
//
vector <Particle *> ReadHVParticles(int &nread){
  nread = 0;
  string cubestartfile="cubestart.dat";
  ifstream incube(cubestartfile.c_str());
  int xstart=-1, ystart=-1,zstart=-1;
  incube>>xstart>>ystart>>zstart;
  PRNT("read_cube",xstart);
  PRNT("read_cube",ystart);
  PRNT("read_cube",zstart);


  if((xstart<0)||(xstart>15)||(xstart<0)||(xstart>15)||(xstart<0)||(xstart>15)){
    cerr<<"[read_cube] Fatal Error: cube start values"<<xstart<<" "<<ystart<<" "<<zstart<<endl;
  }
  vector <Particle*> particles;
  // Reserve max number of particles to avoid memory difficulties
  particles.reserve(MAXSIZE);

  int nxcubes = use_cells;  
  int nycubes = use_cells;
  int nzcubes = use_cells;
  int xfinish = xstart+nxcubes-1;
  int yfinish = ystart+nycubes-1;
  int zfinish = zstart+nzcubes-1;

  int xs,ys,zs;
  xs = xstart;
  ys = ystart;
  zs = zstart;

  cout<<"[read_cube] reading cubes:"<<xs<<" "<<xfinish<<" "<<nxcubes<<" "<<endl;

    for(int i=xs;i<=xfinish;i++) 
      for(int j=ys;j<=yfinish;j++) 
	for(int k=zs;k<=zfinish;k++){
	if(particles.size()>MAXSIZE){
	  cerr<<"[read_cube] Exceeded maximum number of particles"
	      <<particles.size()<<" "<<MAXSIZE<<endl;	  exit(1);
	}
	cout<<"[read_cube] Checking cube"<<i<<" "<<j<<" "<<k<<"...";
	if(ContainsParticles(i,j,k))
	  nread += ReadCube(i,j,k,particles, xstart, ystart, zstart);
	cout<<nread<<endl;
      }
  cout<<"Done."<<endl;
  PRNT("[read_cube] end:",nread);
#endif  
  assert(particles.size()>0);

  return particles;

}
