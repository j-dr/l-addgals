//member functions for vector class
#include <math.h>
#include <iostream.h>
#include "vec3.h"

//default constructor
vec3::vec3() {
  x=0.0;
  y=0.0;
  z=0.0;
  id = 0;
}

//initialize component by component
vec3::vec3(double x1, double x2, double x3) {
  x=x1;
  y=x2;
  z=x3;
}

//set one vec3 equal to another
vec3 vec3::operator = (vec3 vec2) {
  x=vec2.x;
  y=vec2.y;
  z=vec2.z;
  return *this;
}


//return magnitude of vec3
double vec3::magnitude() {
  return sqrt(x*x + y*y + z*z);
}


//add two vec3s
vec3 vec3::operator+ (vec3 vec2) {
  return vec3(x+vec2.x, y+vec2.y, z+vec2.z);
}

//subtract two vec3s
vec3 vec3::operator-(vec3 vec2) {
  return vec3(x-vec2.x, y-vec2.y, z-vec2.z);
}

//multiply a scalar by a vec3
vec3 vec3::operator*(double s) {
  return vec3(s*x, s*y, s*z);
}

//divide a vec3 by a scalar
vec3 vec3::operator/(double s) {
  return vec3(x/s, y/s, z/s);
}

void vec3::print() {
  cout<<"["<<x<<","<<y<<","<<z<<"]";
}

  //should define equality and inequality, but am too lazy at the moment


