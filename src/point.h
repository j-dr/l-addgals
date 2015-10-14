#ifndef point_h
#define point_h
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "simulation.h"
#include "simplefunc.h"
#include "biniostream.h"

using namespace std;
extern Simulation sim;

class Point{
  friend class Particle;
 public:
  Point():x(0),y(0),z(0){};  
  Point(float xx, float yy, float zz):x(xx),y(yy),z(zz){};
  Point(float coord[]):x(coord[0]),y(coord[1]),z(coord[2]){};
    Point(std::vector<float> coord):x(coord[0]),y(coord[1]),z(coord[2]){};
  void Reset(float xx, float yy, float zz){
    x = xx;
    y = yy;
    z = zz;
  };
  void Print()const{cout<<x<<" "<<y<<" "<<z<<endl;};
  void Write(ofstream& file) const{ file<<x<<"\t"<<y<<"\t"<<z<<"\t";}
  void Bin_Write(binfstream& file) const{ file<<x<<y<<z;}
  float X()const{return x;};
  float Y()const{return y;};
  float Z()const{return z;};
  float R()const{return sqrt(x*x+y*y+z*z);};
  float Vector() const{ return sqrt(sqr(x)+sqr(y)+sqr(z));}
  float operator[](int index) const{
    float return_me=0;
    assert( (index >= 0) && (index <3) ); // bounds checking
    if(index==0) return_me = x;
    else if(index==1) return_me = y;
    else if(index==2) return_me = z;
    return return_me;
  }
  float Distance(Point p2) const{
    float dx[3];
    float dist = 0.;
    
    for(int i=0; i<3; i++){
      dx[i] = (*this)[i]-p2[i];
      if (dx[i]<0)
	dx[i] = -dx[i];
      //if (dx[i] > sim.Boxsize()/2){
	//cout << "using periodic boundary conditions"<<dx[i]<<" "<<sim.Boxsize()<<endl;
	//dx[i] = sim.Boxsize() - dx[i];
      //}
      dist += dx[i]*dx[i];    
    }
    dist = sqrt(dist);
    return dist;
  }
  float VerboseDistance(Point p2) const{
    float dx[3];
    float dist = 0.;
    
    for(int i=0; i<3; i++){
      dx[i] = (*this)[i]-p2[i];
      if (dx[i]<0)
	dx[i] = -dx[i];
      if (dx[i] > sim.Boxsize()/2){
	cout << "using periodic boundary conditions"<<dx[i]<<" "<<sim.Boxsize()<<endl;
	dx[i] = sim.Boxsize() - dx[i];
      }
      dist += dx[i]*dx[i];    
      cout<<dx[i]<<" ";
    }
    cout<<dist<<" ";
    dist = sqrt(dist);
    return dist;
  }
 protected:
  float x;
  float y;
  float z;
};

#endif
