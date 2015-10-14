//definitions for a three dimensional vector class
#ifndef vec3_h
#define vec3_h

class vec3 {
 public:
  double x;
  double y;
  double z;
  int id;
public:
  vec3();
  vec3(double, double, double);
  double magnitude();
  vec3 operator=(vec3);
  vec3 operator+(vec3);
  vec3 operator-(vec3);
  vec3 operator-();
  vec3 operator*(double);
  vec3 operator/(double);
  void print();
};

#endif
