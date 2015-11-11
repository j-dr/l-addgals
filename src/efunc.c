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

prefstruct prefs;

void
rng_efunc(const gsl_rng *rng, void *params, double *e1, double *e2)
{
    double e, phi;
    gsl_complex eps;
    eparam *p;
    
    std::cout << "drawing e" << std::endl;
    
    p = (eparam *)params;

    std::cout << "entering do" << std::endl;
    do {
        std::cout << "p->a" << p->a << std::endl;
	std::cout << "p->b" << p->b << std::endl;
	gsl_rng tst = (*rng);
	std::cout << "MSQRT2" << M_SQRT2 << std::endl;
	std::cout << "FUDGE" << FUDGE << std::endl;
	std::cout << "tring different a, b" << std::endl;
	unsigned long int x = gsl_rng_get(rng);
	std::cout << "tring same a, b" << std::endl;
	e = fabs(gsl_ran_exppow(rng, FUDGE * M_SQRT2 * p->a, p->b));
	std::cout << "flat" << std::endl;
	phi = gsl_ran_flat(rng, -M_PI, M_PI);
	std::cout << "dereferenceing e1" << std::endl;
	*e1 = e * cos(2 * phi);
	std::cout << "dereferenceing e2" << std::endl;
	*e2 = e * sin(2 * phi);
	eps = gsl_complex_rect(*e1, *e2);
    } while(gsl_complex_abs(eps) >= prefs.epsmax);
}
