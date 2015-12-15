#include <iostream>
#include "efunc.h"
#include "shapes.h"

/* Because the ellipticity distribution is terminated at |e| < 1, the
   dispersion is a bit too low when drawing it from an analytic model
   that does not fit exactly. Use a fudge factor to get it back up.
   Furthermore, the model was determined on KSB without bias
   correction. Use the Hetterscheidt et al. bias factor of 1/0.85.
 */
#define FUDGE (1.14/0.85)

//prefstruct prefs;

void
rng_efunc(const gsl_rng *rng, void *params, double *e1, double *e2, prefstruct prefs)
{
    double e, phi;
    gsl_complex eps;
    eparam *p;
    
    p = (eparam *)params;
    do {
	gsl_rng tst = (*rng);
	unsigned long int x = gsl_rng_get(rng);
	e = fabs(gsl_ran_exppow(rng, FUDGE * M_SQRT2 * p->a, p->b));
	phi = gsl_ran_flat(rng, -M_PI, M_PI);
	*e1 = e * cos(2 * phi);
	*e2 = e * sin(2 * phi);
	eps = gsl_complex_rect(*e1, *e2);
    } while(gsl_complex_abs(eps) >= prefs.epsmax);
    //std::cout << "e1: " << *e1 << std::endl;
    //std::cout << "e2: " << *e2 << std::endl;
}
