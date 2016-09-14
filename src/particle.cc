#include "particle.h"
#ifdef HEALPIX
#include "healpix_utils.h"
#endif

extern Simulation sim;


void Particle::Print() const{
  cout<<X()<<" "<<Y()<<" "<<Z()<<" "
      <<Xbox()<<" "<<Ybox()<<" "<<Zbox()<<" "
      <<position.x<<" "<<position.y<<" "<<position.z<<" "
      <<Ra()<<" "<<Dec()<<" "<<R()<<" "
      <<Dist8()<<endl;
}



bool Particle::Save() const{

#ifdef FULL_SKY
  if ((position.X()>=-1.0*sim.Boxsize())&&(X()<=sim.Boxsize())&&
      (position.Y()>=-1.0*sim.Boxsize())&&(Y()<=sim.Boxsize())&&
      (position.Z()>=-1.0*sim.Boxsize())&&(Z()<=sim.Boxsize())
      &&(Ra()>=RAMIN)&&(Ra()<=RAMAX)
      &&(Dec()>=DECMIN)&&(Dec()<=DECMAX)
      &&(Zred()>=ZREDMIN)&&(Zred()<=ZREDMAX))
    return true;
  else {
    return false;
  }
#endif

#ifdef BOX_VOLUME
  if ((position.X()>=-0) && (X()<=sim.Boxsize()*BOXFR) &&
      (position.Y()>=-0) && (Y()<=sim.Boxsize()*BOXFR) &&
      (position.Z()>=-0) && (Z()<=sim.Boxsize()*BOXFR))
    return true;
  else {
    return false;
  }
#endif


  float tra = Ra();
  float tr = R();
  float tsindec = Z()/tr;


#ifdef HEALPIX
  long ParticlePixel = 0;
  float phi = tra*PI/180.0;
  float theta = (90.0-Dec())*PI/180.0;
  ang2pix_ring(nSide, theta, phi, &ParticlePixel);

  if ((tra>=RAMIN) && (tra<=RAMAX) &&
      (tsindec>=sinDECMIN) && (tsindec<=sinDECMAX) &&
      (tr>=RMIN_REAL) && (tr<=RMAX_REAL) &&
      (ParticlePixel == PixelNum))
    {
      return true;
    }
  else
    {
      return false;
    }

#endif

  tra = Ra();
  tr = R();
  tsindec = Z()/tr;
  if ((tra>=RAMIN)&&(tra<=RAMAX)
  &&(tsindec>=sinDECMIN)&&(tsindec<=sinDECMAX)
  &&(tr>=RMIN_REAL)&&(tr<=RMAX_REAL))
    return true;
  else {
    return false;
  }
}

bool Particle::Keep() const{
  if (
      (position.X()>0.0)&&(X()<=sim.Boxsize())&&
      (position.Y()>0.0)&&(Y()<=sim.Boxsize())&&
      (position.Z()>0.0)&&(Z()<=sim.Boxsize())&&
      (Ra()>=RAMIN)&&(Ra()<=RAMAX)&&
      (Dec()>=DECMIN)&&(Dec()<=DECMAX)&&
      (Zred()>=ZREDMIN)&&(Zred()<=ZREDMAX))
    return true;
  else {
    cout<<Xbox()<<" "<<Ybox()<<" "<<Zbox()<<" "<<R()<<" "<<Rbox()<<" "<<Zred()<<endl;
    return false;
  }
}
