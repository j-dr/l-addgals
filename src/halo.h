#ifndef halo_h
#define halo_h
#include "cosmology.h"
#ifdef HEALPIX
#include "healpix_utils.h"
#endif

extern Cosmology cosmo;

class Galaxy;

class Halo{
 public:
  Halo(Point pos, Point vel, float m, float zr, int ip, int id, float sig, float rdelta):
    position(pos),velocity(vel),mass(m),zred(zr),np(ip),hid(id),siglos(sig),r200(rdelta){
    ngal=0;
    nmstar=0;
    nbright=0;
    nmid =0;
    ndim=0;
    particle_id = -1;
    host_id = -1;
    dec = angle_const*asin(Z()/R());
    ra = angle_const*atan2(Y(),X());
    if(ra<0) ra = ra+360.;
    centraldistance=rdelta;
    brightest_id=-1;
    central_id=-1;
  };
  void Print()const{cout<<position.X()<<" "<<position.Y()<<" "<<position.Z()<<endl;};
  float M()const{return mass;};
  float PMass()const{return np;};
  float R200()const{
    return r200;};
  float Sig()const{return siglos;};
  float X()const{return position.X();};
  float Y()const{return position.Y();};
  float Z()const{return position.Z();};
  float Vx()const{return velocity.X();};
  float Vy()const{return velocity.Y();};
  float Vz()const{return velocity.Z();};
  Point Position()const{return position;};
  Point Velocity()const{return velocity;};
  int Id()const{return hid;};
  void IncrN(){ngal++;};
  void ResetN(){ngal = 0; nbright=0; nmid=0; ndim=0;};
  int Ngal()const{return ngal;};
  void IncrNbri(){nbright++;};
  int Nbright()const{return nbright;};
  void IncrNmid(){nmid++;};
  int Nmid()const{return nmid;};
  void IncrNdim(){ndim++;};
  int Ndim()const{return ndim;};
  void IncrNms(){nmstar++;};
  int Nmstar()const{return nmstar;};

  double R()const{return sqrt(sqr(X())+sqr(Y())+sqr(Z()));};  
  double Ra()const{return ra;};
  double Dec()const{return dec;};

  double VR()const{return  (Vx()*X()+Vy()*Y()+Vz()*Z())/R();};
#ifdef SNAPSHOT
  double ZredReal()const{return ZREDMIN;};  //z in real space
  double Zred()const{return ZREDMIN;};
#else
  double ZredReal()const{return cosmo.ZofR(R());};  //z in real space
  //with peculiar velocities
  double Zred()const{
    double zz;
    //if(normalization!=SHIFT) 
    zz = cosmo.ZofR(R()); 
    //if(normalization==SHIFT){cout<<"in zred:"<<roffset<<" "<<R()<<" "<<X()<<" "<<Y()<<" "<<Z()<<endl; zz = ZofR(roffset+R())-ZofR(roffset);}
    //if((zz==0.58)&&(simulation==HV)) cout<<"halo:"<<R()<<" "<<X()<<" "<<Y()<<" "<<Z()<<endl;
    return zz+(1+zz)*VR()/cspeed;}; 
#endif


  double Theta(Particle * p)const{
    double pra = p->Ra();
    double pdec = p->Dec();
    double theta;
    //scalar = (size(ra1,/N_Dimen) EQ 0) and (size(ra2,/N_dimen) EQ 0)
    //    double r2d    = 180.0/PI;
    double sindc2 = sin(pdec*d2r);
    double cosdc2 = cos(pdec*d2r);
    double sindc1 = sin(dec*d2r);
    double cosdc1 = cos(dec*d2r);
    double radiff = (ra-pra)*d2r;
    double cosradiff = cos(radiff);
    double cosdist = sindc1*sindc2 + cosdc1*cosdc2*cosradiff;
    assert(cosdist>=-1);
    assert(cosdist<=1);
    //if(cosdist < -1) cosdist = -1.; 
    //if(cosdist > 1) cosdist = 1.; 
    theta = acos(cosdist);
    // this returns theta in radians
    return theta;
  }

  float Dist8()const{return d8;};
  void Dist8(float dist){
    d8 = dist;
  };

  float Mr()const{return mr;};
  void Mr(float mag){
    mr = mag;
  };

  void Central(int gid){
    central_id = gid;
  }

  int Central(){
    return central_id;
  }

  double CentralDistance()const{
    return centraldistance;
  }

  void CentralDistance(float cd){
    centraldistance = cd;
  }

  void Brightest(int gid){
    brightest_id= gid;
  }

  int Brightest()const{
    return brightest_id;
  }
  bool Setgals()const{
    if(brightest_id == -1) return false;
    else return true;
  }

  void Particle(int pid){
    particle_id = pid;
  }

  int Particle(){
    return particle_id;
  }

  void Host(int hid){
    host_id = hid;
  }
  int Host(){
    return host_id;
  }


#ifdef SNAPSHOT
  bool InVol()const{
  if((X()>=0.) && (X()<=sim.Boxsize()) &&
     (Y()>=0.) && (Y()<=sim.Boxsize()) &&
     (Z()>=0.) && (Z()<=sim.Boxsize())) return true;
  else return false;
  }
#else
  bool InVol()const{
#ifdef HEALPIX
    long HaloPixel;
    float phi = Ra()*PI/180.0;
    float theta = (90.0-Dec())*PI/180.0;
    ang2pix_ring(nSide, theta, phi, &HaloPixel);
#endif
    if((Zred()<ZREDMAX)
       &&(Zred()>=ZREDMIN)
       &&(Dec()<DECMAX)
       &&(Dec()>=DECMIN)
       &&(Ra()<RAMAX)
       &&(Ra()>=RAMIN)
#ifdef HVL
       &&(X()<use_cells*187.5)
       &&(Y()<use_cells*187.5)
       &&(Z()<use_cells*187.5)
#endif
#ifdef HEALPIX
       &&(HaloPixel==PixelNum)
#endif
       )
      return true;
    else return false;
  }
#endif

 private:
  Point position;  // halo position (3-d point)
  Point velocity;  // halo velocity (3-d point)
  float mass;
  float zred;
  float np;
  int hid;
  int ngal;
  int nmstar;
  int nbright;
  int nmid;
  int ndim;
  float siglos;
  float r200;
  double ra;
  double dec;
  float d8;
  float mr;
  float centraldistance;
  int central_id;
  int brightest_id;
  int particle_id;
  int host_id;
};

#endif
