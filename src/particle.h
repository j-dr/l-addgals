#ifndef particle_h
#define particle_h
#include <cmath>        //for trig
#include <iostream>     //for cout
#include "constants.h"
#include "box.h"
//#include "simplefunc.h" //for sqr
#include "pi.h" //for pi
#include "point.h"
#include "cosmo.h"
#include "functions.h"
#include "cosmology.h"
#include "global_vars.h"

using namespace std;
extern Cosmology cosmo;

struct hrow{
  long hid;
  double rhalo;
  double mvir;
  double rvir;
};


class Particle{
  friend class Galaxy;
 public:
  Particle(){gid=-1; hid=-1; pid=-1; mhost=0.; rvir=-1; rhalo=-1.;};
  Particle(Point pos, Point vel, float dens):
  position(pos),velocity(vel),density(dens){ gid=-1; hid=-1; pid=-1; mhost=0.; rvir=-1.; rhalo=-1.;};

  void Hid(int id){hid=id;};
  int Hid()const{return hid;};
  void Mhost(float mass){mhost=mass;};
  float Mhost()const{return mhost;};
#ifdef LONG64_IDS
  void Pid(long int id){pid=id;};
  int Pid()const{return pid;};
#else
  void Pid(int id){pid=id;};
  long int Pid()const{return pid;};
#endif
  void PosPrint()const{cout<<X()<<" "<<Y()<<" "<<Z()<<endl;};
  void PosAssign(Point tpos){position=tpos;};
  void VelAssign(Point tvel){velocity=tvel;};
  void DensAssign(float tdens){density=tdens;};
  void HaloAssign(hrow hdata){
      hid = hdata.hid;
      rhalo = hdata.rhalo;
      rvir = hdata.rvir;
      mhost = hdata.mvir;
  }
  void RVir(float rv){rvir=rv;};
  void RHalo(float rh){rhalo=rh;};
  void MVir(float mv){mhost=mv;};
  void VelPrint()const{cout<<Vx()<<" "<<Vy()<<" "<<Vz()<<endl;};
  void Print() const;
  void Write(ofstream &file)const{
    file<<X()<<" "<<Y()<<" "<<Z()<<" "<<Vx()<<" "<<Vy()<<" "<<Vz()<<" "<<" "<<Zred()<<" "<<Den()<<" "<<Pid()<<endl;
  }
  void Pos2Write(ofstream &file)const{
    file<<X()<<" "<<Y()<<" "<<Z()<<" "<<Ra()<<" "<<Dec()<<" "<<Zred()<<endl;
  }

  bool Save() const;
  bool Keep() const;
  Point Position() const{Point pos(X(),Y(),Z()); return pos;};
  double Distance(Point p2)const{    return Position().Distance(p2);};
  double Distance(Particle *p2)const{    return Position().Distance(p2->Position());};
  //  double Dist2D(Point p2)const{    return Position().Dist2D(p2);};
  //double Dist2D(Particle *p2)const{    return Position().Dist2D(p2->Position());};
  float Dist8()const{return density;};
  Point BoxPosition(){return position;};
  Point Velocity(){return velocity;};

  double X()const{return sim.LengthUnit()*Xbox();};
  double Y()const{return sim.LengthUnit()*Ybox();};
  double Z()const{return sim.LengthUnit()*Zbox();};
  double Xbox()const{return position.x;};
  double Ybox()const{return position.y;};
  double Zbox()const{return position.z;};

  float Vx()const{return 1.0*velocity.x;};
  float Vy()const{return 1.0*velocity.y;};
  float Vz()const{return 1.0*velocity.z;};
  float Den()const{return density;};
  double Rbox()const{return sqrt(sqr(position.x-0.5)+sqr(position.y-0.5)+sqr(position.z-0.5));};
  double R()const{return sqrt(sqr(X())+sqr(Y())+sqr(Z()));};

  double Ra()const{double ra = angle_const*atan2(Ybox(),Xbox()); if(ra<0) return ra+360.; else return ra;}
  double Dec()const{return angle_const*asin(Z()/R());};
  double VR()const{return  (Vx()*X()+Vy()*Y()+Vz()*Z())/R();};
  float RHalo()const{return rhalo;};
  float RVir()const{return rvir;};
  float MVir()const{return mhost;};
  double ZredReal()const{
#ifdef SNAPSHOT
  return sim_redshift;
#else
  float zz = cosmo.ZofR(R());
  //if(SHIFT) zz = ZofR(roffset+R())-ZofR(roffset);
  return zz;
#endif
  };
  //with peculiar velocities
  double Zred()const{
#ifdef SNAPSHOT
  return sim_redshift;
#else
  double zz = cosmo.ZofR(R());
  return zz+(1+zz)*VR()/cspeed;
#endif
  };
  //this is the angle from the line of sight
  //this is from the idl program (my)gcirc

  void MakeGal(int id){
    gid=id;
  }

  void SetZred(float zGal){
    double rGal = cosmo.RofZ(zGal);
    double rfac = rGal/R();
    position.x *= rfac;
    position.y *= rfac;
    position.z *= rfac;
  }

  void UnMakeGal(){
    gid=-1;
  }

  bool IsGal()const{
    if(gid>=0) return true;
    return false;
  }

  int Gid()const{
    return gid;
    }
  bool InHalo()const{if(hid==-1) return false; else return true;};//return inhalo;};
  void SetHaloId(int id){
    hid=id;
  }
  void DeleteOutofRangeParticle(Particle*& p)
  {
    if( (p->Xbox()>=XMIN)&&(p->Xbox()<=XMAX)&&
	(p->Ybox()>=YMIN)&&(p->Ybox()<=YMAX)&&
	(p->Zbox()>=ZMIN)&&(p->Zbox()<=ZMAX)){

      if((p->Ra()>=RAMIN)&&(p->Ra()<=RAMAX)&&
	 (p->Dec()>=DECMIN)&&(p->Dec()<=DECMAX)&&
	 (p->Zred()>=ZREDMIN)&&(p->Zred()<=ZREDMAX)){
	//cout<<p->X()<<" "<<p->Y()<<" "<<p->Z()<<" "<<p->Ra()<<" "<<p->Dec()<<" "<<p->Zred()<<endl;
      }
    }
    else{
      cout<<"[hv:DeleteOutofRangeParticle] deleting: "<<p->Xbox()<<" "<<p->X()<<" "<<p->Y()<<" "<<p->Z()<<" "<<p->Ra()<<" "<<p->Dec()<<" "<<p->Zred()<<endl;
      delete p;
      p=0;
    }
  }

 private:
  Point position;  // particle position (3-d point)
  Point velocity;  // particle velocity (3-d point)
  float density;  //eigth nearest particle
  int gid;
  int hid;
  float mhost;
  float rvir;
  float rhalo;
#ifdef LONG64_IDS
  long int pid;
#else
  int pid;
#endif
};




#endif
