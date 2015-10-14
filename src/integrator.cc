// integrator.cc
//
// Implementation of numerical integrator routines
// Based on fortran code my Mike Gross
// Uses fifth order Runge-Kutta with adaptive step size
//
// Patrik Jonsson
//

#include "./integrator.h"
#include <cmath>
#include <iostream>
#include <stdarg.h>

using namespace std;
// Constructor for an integrator
// An integrator has a specified function, accuracy and initial step size
integrator::integrator(double(*i)(double,double*), int n, double a, double s)
{
  integrand = i;
  parameters = 0;
  accuracy = a;
  startstep = s;
  npar = n;

  parameters = new double[npar]; // Allocate array of parameters
}

// Call to Integrate returns integral of function defined in constructor
// between a and b. The ellipsis are n additional parameters propagated to
// the integrand and had better conform to its syntax.
double integrator::Integrate(double a, double b, ...)
{
  const int maxsteps = 100000000;

  double dx, dxnext;
  int Nstep;

  if(npar>0) {
    // Find parameters given to integrand and put them in the array
    va_list ap;
    va_start(ap,b);
    for(int i=0;i<npar;i++)
      parameters[i] = va_arg(ap,double);
    va_end(ap);
  }

  // Initialise members
  x = a;
  y = 0;
  dx = startstep;
  Nstep = 0;

  while( ( ((x-b)*(b-a))<0 ) && (Nstep<maxsteps)) {
    Nstep++;
    dydx = integrand(x,parameters);

    // yscale is used to monitor accuracy. This general purpose choice
    // can be modified if need be.

    yscale = ( (fabs(y)+fabs(dx*dydx))>1e-16 ? fabs(y)+fabs(dx*dydx) : 1e-16 );

    if( ((x+dx-b)*(x+dx-a)) > 0) // Stepsize overshoots, decrease it
      dx = b-x;

    dx = Makestep(dx);
  }

  if(Nstep==maxsteps)
    cout << "integrator::integrate: Warning, failed to converge\n";

  return y;
}

// 5th order Runge-Kutta step
double integrator::Makestep(double htry)
{
  const double safety=0.9, pgrow=-0.2, pshrink=-0.25, errcon=1.89e-4;

  double errmax, h, yerr, ytemp;
  int underflow;

  h = htry; // Set stepsize to initial accuracy
  errmax = 10;
  underflow = 0;

  while( (errmax>1) && (!underflow) ) {
    ytemp = Runge5(h, yerr); // Take a step

    errmax = fabs(yerr/yscale)/accuracy; // Scale relative to required accuracy

    if(errmax>1) { // Truncation error too large, reduce h
      double htemp = safety*h*pow(errmax,pshrink);
      h = sign( (fabs(htemp)>0.1*fabs(h)) ? fabs(htemp) : 0.1*fabs(h) , h);
      double xnew = x+h;

      if( (xnew==x) && (!underflow) ) {
	underflow = 1;
	cout << "integrator::Makestep: Warning, stepsize underflow\n";
      }
    }
  }


  x+=h;
  y=ytemp;

  // Step succeeded. Estimate and return size of next step

  if(errmax>errcon)
    return safety*h*pow(errmax,pgrow);
  else
    return 5.0*h; // Not more than a factor of 5 increase
}


// Advance solution using 5th order C-K R-K
// Returns y-value and error in reference argument yerr
double integrator::Runge5(double h, double& yerr)
{
  const double a3=0.3, a4=0.6, a5=1.0, a6=0.875;
  const double c1=37.0/378, c3=250.0/621, c4=125.0/594, c6=512.0/1771;
  const double dc1= c1-2825.0/27648, dc3= c3-18575.0/48384;
  const double dc4= c4-13525.0/55296, dc5=-277.0/14336;
  const double dc6= c6-0.25;

  double ak3,ak4,ak5,ak6;

  ak3 = integrand(x+a3*h, parameters);
  ak4 = integrand(x+a4*h, parameters);
  ak5 = integrand(x+a5*h, parameters);
  ak6 = integrand(x+a6*h, parameters);

  // Estimate error as difference between fourth and fifth order
  yerr = h*(dc1*dydx + dc3*ak3 + dc4*ak4 + dc5*ak5 + dc6*ak6);

  // Estimate fifth-order value
  return (y + h*(c1*dydx + c3*ak3 + c4*ak4 + c6*ak6) );


}

// f77 function sign
double integrator::sign(double val, double sgn)
{
  double s = sgn/fabs(sgn);
  return fabs(val)*s;
}

    
    
  

