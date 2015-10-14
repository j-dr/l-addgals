#ifndef clusters_h
#define clusters_h
#include "cosmology.h"

extern Cosmology cosmo;

class Cluster{
 public:
  Cluster(Point pos, Point vel, int n, float rra, float ddec, float zz, float r, int id):
    position(pos),velocity(vel),ngals(n), ra(rra), dec(ddec), zred(zz), r200(r), hid(id){
  };
  void Print()const{cout<<position.X()<<" "<<position.Y()<<" "<<position.Z()<<endl;};
  float R200()const{
    return r200;};
  //  float Sig()const{return siglos;};
  float X()const{return position.X();};
  float Y()const{return position.Y();};
  float Z()const{return position.Z();};
  Point Position()const{return position;};
  float Vx()const{return velocity.X();};
  float Vy()const{return velocity.Y();};
  float Vz()const{return velocity.Z();};
  Point Velocity()const{return velocity;};
  int Id()const{return hid;};
  //  int Ngal()const{return ngal;};
  int Ngals()const{return ngals;};
  double R()const{return sqrt(sqr(X())+sqr(Y())+sqr(Z()));};  
  double Ra()const{return ra;};
  double Dec()const{return dec;};
  double VR()const{return  (Vx()*X()+Vy()*Y()+Vz()*Z())/R();};
  double ZredReal()const{return cosmo.ZofR(R());};  //z in real space
  //with peculiar velocities
  double Zred()const{
    return zred;
  }
  double Theta(Particle * p)const{
    double pra = p->Ra();
    double pdec = p->Dec();
    double theta;
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
 private:
  Point position;  // halo position (3-d point)
  Point velocity;  // halo velocity (3-d point)
  int ngals;
  //float siglos;
  double ra;
  double dec;
  double zred;
  float r200;
  int hid;
  };

#endif
