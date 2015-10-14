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
// boolian function to determine whether cell ijk contains any particles
// within requested limits
//
bool ContainsParticles(int i, int j, int k){
  float hb = BOXES/2.0;
  float rmin = sqrt(sqr(i-hb)+sqr(j-hb)+sqr(k-hb))/hb;
  float rmax = sqrt(sqr(i+1-hb)+sqr(j+1-hb)+sqr(k+1-hb))/hb;
  cout<<"[read_cube]"<<rmin<<" "<<rmax<<endl;
  PRNT("read_cube/Contains",1.0*i/BOXES);
  PRNT("read_cube/Contains",1.0*j/BOXES);
  PRNT("read_cube/Contains",1.0*k/BOXES);
  if((((float)i/BOXES)>=XMIN)&&((i+1.0)/BOXES<=XMAX)&&
     (((float)j/BOXES)>=YMIN)&&((j+1.0)/BOXES<=YMAX)&&
     (((float)k/BOXES)>=ZMIN)&&((k+1.0)/BOXES<=ZMAX))//&&
     // (rmin>=RMIN)&&(rmax<=RMAX))
    return true;
  else{
    cout<<((float)i/BOXES)<<" "<<XMIN<<" ";
    cout<<(i+1.0)/BOXES<<" "<<XMAX<<endl;
   cout<<"[read_cube] No particles in cube"<<i<<j<<k<<"; not reading"<<endl;
    return false;
    }
}


//
// boolian function to determine whether (x,y,z) is in the volume requested
//

//bool PointInVol(Point p){
  //***written for cube centered at origin going positive from 8***//
  //if((p.X()>(xstart-8)*LengthUnit/BOXES)&&(p.X()<(xstart-8)+use_cells*LengthUnit/BOXES)&&
// (p.Y()>(ystart-8)*LengthUnit/BOXES)&&(p.Y()<(ystart-8)+use_cells*LengthUnit/BOXES)&&
//   (p.Z()>(zstart-8)*LengthUnit/BOXES)&&(p.Z()<(zstart-8)+use_cells*LengthUnit/BOXES))
//  return true;
//else return false;
//}


//
// boolian function to determine whether halo is in the (ra,dec,z) volume requested
//
/*bool HaloInVol(Halo * h){
  if((h->Zred()<ZREDMAX)
     &&(h->Dec()<DECMAX)
     &&(h->Dec()>DECMIN)
     &&(h->Ra()>RAMIN)
     &&(h->Ra()<RAMAX)
     &&(h->X()<use_cells*187.5)
     &&(h->Y()<use_cells*187.5)
     &&(h->Z()<use_cells*187.5))
    return true;
  else return false;
}
*/

//
// Converts short int position into box length units
//
float PositionConvert(short int x, int cube, int xs){
  float pos;
  if(simulation==HV){
    pos = (((float) x)/((float) MMAXINT)+((float) cube))/((float) BOXES);
  //PRNTV(x);
  //PRNTV(pos);
  //PRNTV((((float) x)/((float) MMAXINT)+((float) cube))/((float) BOXES));
  //    cout<<x<<" "<<pos<<" "<<p2<<endl;
  //<<endl;
  //PRNTV(xs);
  
    if(normalization == SHIFT){
      float offset_cells = xs-8;
      if(xs<=7) offset_cells = 7-xs;
      pos = pos-0.0625*offset_cells;
    }
    if(box!=HVPOW){
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
  }
  else{
    //other non-hv sims.  where does this apply?
    pos = (((float) x)/((float) MMAXINT)+((float) cube))/((float) BOXES);
    //pos = x;
    cout<<x<<" ";
    //pos = ((float) x + (float (MMAXINT)))/(2*(float (MMAXINT)));
    //pos = (pos+((float) cube))/((float) BOXES);
    cout<<pos<<endl;
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
  cout<<"[read_cube] looking for cube "<<i<<" "<<j<<" "<<k<<endl; 
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
  else cerr<<"reading file '"<<fstrp<<"'"<<endl;
  if (vfile.fail()) {
    cerr<<"error: cannot open file '" <<fstrv<<"'"<<endl;
    return 0;
  }
  else cout<<"reading file '"<<fstrv<<"'"<<endl;
  
  if (hfile.fail()) {
    cerr<<"error: cannot open file '" <<fstrh<<"'"<<endl;
    return 0;
  }
  else cout<<"reading file '"<<fstrh<<"'"<<endl;

  if (hidfile.fail()) {
    cerr<<"error: cannot open file '" <<fstrhid<<"'"<<endl;
    return 0;
  }
  else cout<<"reading file '"<<fstrhid<<"'"<<endl;

  cout<<"[read_cube] reading files "<<fstrp<<" "<<fstrv<<" "<<fstrh<<" "<<fstrhid<<endl;

  //#ifndef SPARC
  //convert endian-ness on intel & dec machines
  //you need this on ghostwheel.
  //ghostwheel is little endian.  the files are big endian.
  //therefore you need to read them as such.

#ifdef LITTLEENDIAN
  pfile.bigendian();
  vfile.bigendian();
  hfile.bigendian();
#endif


  //#endif

  short int tmp_x, tmp_y, tmp_z, tmp_dens;
  unsigned int np = 0;
  unsigned int nr = 0;

  int hid;

  /*
  if(normalization == SHIFT){
    cout<<"[read_cube] ***  ATTENTION  ***"<<endl;
    cout<<"[read_cube] Converting to a less evolved universe!"<<endl;
    assert(xstart > 8);
    assert(ystart==xstart);
    assert(zstart==xstart);
    float offset = xstart-8;
  //  cout<<"[read_cube] Offset by "<<ZofR(offset*187.5)<<endl;
  }
  */
  //char fstrhout[20];
  //sprintf(fstrhout,"/data2/risa/tmp_hvdens/h.%02i.%02i.%02i.txt",i,j,k);
  //ofstream hout(fstrhout);

  while((pfile)&&(vfile)&&(hfile)){
  //while((pfile)&&(vfile)){
    //** Read 3 short ints from particle file
    //** Convert to box positions as they are passed to constructor
    pfile>>tmp_x>>tmp_y>>tmp_z;
    Point xx(PositionConvert(tmp_x,i,xstart), 
	     PositionConvert(tmp_y,j,ystart), 
	     PositionConvert(tmp_z,k,zstart));
    //cout<<tmp_x<<" "<<tmp_y<<" "<<tmp_z<<" "<<xx.X()<<" "<<xx.Y()<<" "<<xx.Z()<<endl;
    //if((xx.X()<0.5)||(xx.X()>1))     xx.Print();
    //if((xx.Y()<0.5)||(xx.Y()>1))     xx.Print();
    //if((xx.Z()<0.5)||(xx.Z()>1))     xx.Print();

    //** Read 3 short ints from velocity file, save as ShortPoint
    vfile>>tmp_x>>tmp_y>>tmp_z;
    //ShortPoint vv(tmp_x, tmp_y, tmp_z);
    Point vv(tmp_x, tmp_y, tmp_z);

    //** Read 1 short int from density file
    hfile>>tmp_dens;
    //cout<<tmp_dens<<endl;
    //ddfile<<tmp_dens<<" ";
    nr++;

    float hminabs = 1. / (16.*2.*MMAXINT);
    float sfaci =  (log(1.) - log(hminabs))/  (2.*MMAXINT);
    float dist = hminabs * exp( sfaci* (tmp_dens +MMAXINT) );
    //    cout<<sim.LengthUnit()<<" ";
    float dist8= dist*sim.LengthUnit();
    hidfile>>hid;   


#ifdef CATALOG
    //    cout<<"Catalog"<<endl;
    Particle * particle = new Particle(xx, vv, dist8);
    if(hid>=0){
      particle->Hid(hid);

    }
    else
      cout<<np<<" "<<hid<<endl;
    //      if((particle->Xbox()<0) && (xstart > 7))
    //if((particle->Xbox()<0))
    //	cout<<"[read_cube] X:"<<particle->Xbox()<<" "<<particle->X()<<endl;
    if(particle->Save()){
      particles.push_back(particle);
      np++;
      if(particle->Zred()>zmax->GetVal())
	zmax->SetVal(particle->Zred());
    }
    else{
      //cout<<"p outside"<<particle->Ra()<<" "<<particle->Dec()<<" "
      //  <<particle->X()<<" "<<particle->Y()<<" "<<particle->Z()<<endl;
      delete particle;
    }

    //    else
    //cout<<"negative hid"<<hid<<endl;
#endif
#ifdef HALOASSIGN
    //    cout<<"Halo assign"<<endl;
    Particle * particle = new Particle(xx, vv, dist8);
    particles.push_back(particle);
    //	cout<<particle->Xbox()<<" "<<particle->Ybox()<<" "<<particle->Zbox()<<" ";
    //cout<<particle->X()<<" "<<particle->Y()<<" "<<particle->Z()<<endl;
    np++;
#endif


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
#ifdef HVL
  if(box==HVLIGHT){
    if(xstart <8) {
      xfinish = xstart;
      xs = xfinish-nxcubes+1;
    } else xs=xstart;
    
    if(ystart <8) {
      yfinish = ystart;
      ys = yfinish-nycubes+1;
    } else ys = ystart;
    
    if(zstart <8) {
      zfinish = zstart;
      zs = zfinish-nzcubes+1;
    } else zs = zstart;
  }
  else{
    xs = xstart;
    ys = ystart;
    zs = zstart;
  }

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
