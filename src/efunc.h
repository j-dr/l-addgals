#ifndef _EFUNC_H_
#define _EFUNC_H_ 1

#include <math.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>
#include "shapeconfig.h"

/*!
 * \struct eparam
 * \brief PDF parameters for ellipticity distribution
 */
struct eparam {
    double a, b;
};
typedef struct eparam eparam;

double efunc(double e, void *);
double efunc_norm(double, double);
void rng_efunc(const gsl_rng *, void *, double *, double *, prefstruct);

#endif
