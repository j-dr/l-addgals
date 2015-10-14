#include "myrand.h"

//returns a random integer between 1 and m
int randint(int m){
  assert (m>0);
  int return_me =  static_cast <int> (ceil(drand48()*m));
  //should be unnecessary, but just as a sanity check...
  assert ((return_me>=1));
  assert ((return_me<=m));
  return return_me;
}

//returns a random integer between m1 and m2, inclusive
int randint(int m1, int m2){
  assert(m1>=0);
  assert(m2>=m1);
  int ri = randint(m2-m1+1);
  return ri+m1-1;
}

//converts a float/double to an integer
// floor with prob: f-i.
int randround(double ff){
  double rn = drand48();
  int return_int = static_cast <int> (floor(ff));
  double remain = ff-static_cast <int> (floor(ff));
  if(rn < remain) return_int++;
  return return_int;
}

//returns true with probability p
bool randbool(double p){
  double rn = drand48();
  if(rn<p) return true;
  else return false;
}
