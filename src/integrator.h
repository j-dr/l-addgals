// integrator.h
//
// Declaration of numerical integrator routines
// Based on fortran code my Mike Gross
// Uses fifth order Runge-Kutta with adaptive step size
//
// Patrik Jonsson
//

/*
One integrator object is initialised to integrate a specific function
to a specified accuracy with a specified starting step.
The syntax of the integrating function is
double integrand(double x, double* parameters)
where the parameters are an arbitrary number of parameters that are not
integrated over.

Integration is performed with
integral = i.Integrate(lower, upper, ... )
where the ellipsis should be the same number of parameters that the integrand
is taking.
*/

#ifndef __integrator__
#define __integrator__

class integrator {
private:
  double (*integrand)(double,double*); // Pointer to integrand function
  int npar;                            // Number of integrand parameters
  double* parameters;                  // Integrand parameters
  double accuracy;                     // Required accuracy;
  double startstep;                    // First-guess step size

  double x,y,dydx,yscale;              // Integration variables

  double Runge5(double h,double& yerr);
  double Makestep(double htry);
  double sign(double val, double sgn);

public:
  integrator( double(*i)(double, double*), int n, double a, double s);

  double Integrate(double a, double b, ...); // Make integration

};

#endif

