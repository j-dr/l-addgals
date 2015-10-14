#ifndef myrand_h
#define myrand_h
#include <cassert>
#include <cmath>
#include <cstdlib>    //for drand48()

int randint(int m);
int randint(int m1, int m2);
int randround(double dd);
bool randbool(double p);

#endif
