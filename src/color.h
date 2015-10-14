#ifndef color_h
#define color_h
#include <cassert>    
#include <iostream>   //for i/o
#include <fstream>    //for file i/o
#include <vector>     //for vectors
#include <cstdlib>    //for exit()
#include <cmath>      //for floor()
#include <string>     //for strings
#include "myrand.h"    //for randint
#include "constants.h" //for Rmin
#include "cosmo.h"    //for app/abs conversion



//** Choose colors from real galaxies, based on either
//** luminosity or local density 
//#define COLORDENS
//#define COLORLUM
//#define COLORDLUM
// old options removed

using std::cout;

#define SEDDLUM

class GalSED{
 public:
 GalSED(double absr, double dd, int _id, int _catid):Mr(absr), dens(dd), id(_id), catid(_catid), red(0){};
 GalSED(double absr, double dd, int _id, int _catid, int _red):Mr(absr), dens(dd), id(_id), catid(_catid), red(_red){};
  double MR()const{return Mr;};
  double Dens()const{return dens;};
  int Id()const{return id;};
  int CatId()const{return catid;};
  int Red()const{return red;};
 private:
  double Mr;
  double dens;
  int id;
  int catid;
  int red;
};


#ifdef DR3
class GalSED{
 public:
  GalSED(double absr, double dd):Mr(absr), dens(dd){};
  double MR()const{return Mr;};
  double Dens()const{return dens;};
  int Id()const{return id;};
  void Id(int _id){id=_id;};
 private:
  //vector <double> sed;
  double Mr;
  double dens;
  int id;
};                         
#endif

#ifdef DR2

class SEDTuple{
public:
  SEDTuple(){};
  SEDTuple(float tc1, float tc2, float tc3, float tMr):
    c1(tc1),c2(tc2),c3(tc3),mr(tMr){};
  void Print(){cout<<c1<<" "<<c2<<" "<<c3<<endl;};
  void Write(ofstream &file){file<<c1<<" "<<c2<<" "<<c3<<" ";};
  float operator[](int index) const{
    float return_me=0;
    assert( (index >= 0) && (index <4) ); // bounds checking
    if(index==0) return_me = c1;
    else if(index==1) return_me = c2;
    else if(index==2) return_me = c3;
    else if(index==3) return_me = mr;
    return return_me;
  }
  //const float Mr(){return mr;);
  
private:
  float c1;
  float c2;
  float c3;
  float c4;
  float mr;
};

#endif

#ifdef DR3
class SEDTuple{
 public:
  SEDTuple(){};
  SEDTuple(float tc1, float tc2, float tc3, float tMr,  vector <int> sid):
  c1(tc1),c2(tc2),c3(tc3),mr(tMr), id(sid){};
  void Print(){cout<<mr<<" "<<c1<<" "<<c2<<" "<<c3<<" "<<endl;};
  void Write(ofstream &file){file<<c1<<" "<<c2<<" "<<c3<<" ";};
  float operator[](int index) const{
    float return_me=0;
    assert( (index >= 0) && (index <4) ); // bounds checking
    if(index==0) return_me = c1;
    else if(index==1) return_me = c2;
    else if(index==2) return_me = c3;
    else if(index==3) return_me = mr;
    return return_me;
  }
  //const float Mr(){return mr;);
  
 private:
  float c1;
  float c2;
  float c3;
  float mr;
  vector <int> id;
};
#endif

/*
class SEDFiveTuple{
public:
  SEDFiveTuple(){};
  SEDFiveTuple(float tc1, float tc2, float tc3, float tc4, float tMr):
    c1(tc1),c2(tc2),c3(tc3),c4(tc4),mr(tMr){};
  void Print(){cout<<c1<<" "<<c2<<" "<<c3<<" "<<c4<<endl;};
  void Write(ofstream &file){file<<c1<<" "<<c2<<" "<<c3<<" "<<c4<<" ";};
  const float operator[](int index) const{
    float return_me=0;
    assert( (index >= 0) && (index <5) ); // bounds checking
    if(index==0) return_me = c1;
    else if(index==1) return_me = c2;
    else if(index==2) return_me = c3;
    else if(index==3) return_me = c4;
    else if(index==4) return_me = mr;
    return return_me;
  }
  //const float Mr(){return mr;);
  
private:
  float c1;
  float c2;
  float c3;
  float c4;
  float mr;
};
 */
#endif
