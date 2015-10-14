#include <math.h>

/*************************    Program 2.13   ****************************/
/*                                                                      */
/************************************************************************/
/* Please Note:                                                         */
/*                                                                      */
/* (1) This computer program is written by Tao Pang in conjunction with */
/*     his book, "An Introduction to Computational Physics," published  */
/*     by Cambridge University Press in 1997.                           */
/*                                                                      */
/* (2) No warranties, express or implied, are made for this program.    */
/*                                                                      */
/************************************************************************/
int iseed = 1;

double ranf()
/* Uniform random number generator x(n+1)= a*x(n) mod c
   with a = pow(7,5) and c = pow(2,31)-1.
   Copyright (c) Tao Pang 1997. */
{
const int ia=16807,ic=2147483647,iq=127773,ir=2836;
int il,ih,it;
double rc;
//extern int iseed;
ih = iseed/iq;
il = iseed%iq;
it = ia*il-ir*ih;
if (it > 0)
  {
  iseed = it;
  }
else
  {
iseed = ic+it;
  }
rc = ic;
return iseed/rc;
}

void grnf (double *x,double *y)
/* Two Gaussian random numbers generated from two uniform
   random numbers.  Copyright (c) Tao Pang 1997. */
//double *x,*y;
{
double pi,r1,r2;
double ranf();

pi =  4*atan(1);
r1 = -log(1-ranf());
r2 =  2*pi*ranf();
r1 =  sqrt(2*r1);
*x  = r1*cos(r2);
*y  = r1*sin(r2);
}

double normal_random(float mean, float stddev) {
  double x,y;
  grnf(&x, &y);
  return (x*stddev + mean);
}
