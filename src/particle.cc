#include "particle.h"
#ifdef HEALPIX
#include "chealpix.h"
#endif

extern Simulation sim;

/*
float Point::Distance(Point p2) const{
  float dx[3];
  float dist = 0.;

  for(int i=0; i<3; i++){
    dx[i] = (*this)[i]-p2[i];
    if (dx[i]<0)
      dx[i] = -dx[i];
    //if (dx[i] > BOXSIZE/2){
      //cout << "using periodic boundary conditions\n";
    // dx[i] = BOXSIZE - dx[i];
    //}
    //  cout<<dx[i]<<endl;
    dist += dx[i]*dx[i];
  }
  dist = sqrt(dist);
  // cout<<dist<<endl;
  return dist;
}
*/
/*
float Point::Dist2D(Point p2) const{
  float dx[2];
  float dist = 0.;

  for(int i=0; i<2; i++){
    dx[i] = (*this)[i]-p2[i];
    if (dx[i]<0)
      dx[i] = -dx[i];
    dist += dx[i]*dx[i];
  }
  dist = sqrt(dist);
  return dist;
}

*/
/*
float Point::Dist2D(Point p2) const{
  return sqrt(sqr(x-p2.X())+sqr(y-p2.Y()));
}
*/

//float Point::ProjDist(Point p2) const{
// return sqrt(sqr(x-p2.X())+sqr(y-p2.Y()));
//}

void Particle::Print() const{
  cout<<X()<<" "<<Y()<<" "<<Z()<<" "
      <<Xbox()<<" "<<Ybox()<<" "<<Zbox()<<" "
      <<position.x<<" "<<position.y<<" "<<position.z<<" "
    //  <<Vx()<<" "<<Vy()<<" "<<Vz()<<" "
    //      <<Zred()<<" "
      <<Ra()<<" "<<Dec()<<" "<<R()<<" "
      <<Dist8()<<endl;
}



bool Particle::Save() const{
  //return true;
  //cout<<Xbox()<<" "<<Ybox()<<" "<<Zbox()<<" ";
  //cout<<Ra()<<" "<<Dec()<<" "<<Zred()<<endl;
#ifdef FULL_SKY
  if ((position.X()>=-1.0*sim.Boxsize())&&(X()<=sim.Boxsize())&&
      (position.Y()>=-1.0*sim.Boxsize())&&(Y()<=sim.Boxsize())&&
      (position.Z()>=-1.0*sim.Boxsize())&&(Z()<=sim.Boxsize())
      &&(Ra()>=RAMIN)&&(Ra()<=RAMAX)
      &&(Dec()>=DECMIN)&&(Dec()<=DECMAX)
      &&(Zred()>=ZREDMIN)&&(Zred()<=ZREDMAX))
    return true;
  else {
    //        cout<<"bad p"<<sim.Boxsize()<<" "
    //        cout
    //<<sim.LengthUnit()<<" "
    //<<position.X()<<" "
    //<<position.Y()<<" "
    //<<position.Z()<<" "
    //<<X()<<" "<<Y()<<" "<<Z()<<" "<<endl;
    // <<R()<<" "<<Rbox()<<" "<<ZredReal()<<endl;
    //<<Ra()<<" "<<Dec()<<" "<<Zred()<<endl;
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
  //float theta = (90.0-asin(tsindec))*PI/180.0;
  float theta = (90.0-Dec())*PI/180.0;
  ang2pix_ring(nSide, theta, phi, &ParticlePixel);
  //cout<<ParticlePixel<<" "<<PixelNum<<" "<<RAMIN<<" "<<tra<<" "<<RAMAX<<" "<<sinDECMIN<<" "<<tsindec<<" "<<sinDECMAX<<" "<<RMIN_REAL<<" "<<tr<<" "<<RMAX_REAL<<endl;
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
  //if ((position.X()>=0.0)&&(X()<=sim.Boxsize())&&
  //(position.Y()>=0.0)&&(Y()<=sim.Boxsize())&&
  //(position.Z()>=0.0)&&(Z()<=sim.Boxsize())
  //&&(Ra()>=RAMIN)&&(Ra()<=RAMAX)
  //&&(Dec()>=DECMIN)&&(Dec()<=DECMAX)
  //&&(Zred()>=ZREDMIN)&&(Zred()<=ZREDMAX))
  //return true;
  else {
    //        cout<<"bad p"<<sim.Boxsize()<<" "
    //        cout
    //<<sim.LengthUnit()<<" "
    //<<position.X()<<" "
    //<<position.Y()<<" "
    //<<position.Z()<<" "
    //<<X()<<" "<<Y()<<" "<<Z()<<" "<<endl;
    // <<R()<<" "<<Rbox()<<" "<<ZredReal()<<endl;
    //<<Ra()<<" "<<Dec()<<" "<<Zred()<<endl;
    return false;
  }
}

bool Particle::Keep() const{
  if (//(position.x>=XMIN)&&(position.x<=XMAX)&&
      //(position.y>=YMIN)&&(position.y<=YMAX)&&
      //(position.z>=ZMIN)&&(position.z<=ZMAX)&&
      //(R()>=RMIN)&&(R()<=RMAX)&&
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

