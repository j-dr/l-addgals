/*!\file gno.c
 * \brief Generalized Normal Distribution
 */
#include <iostream>
#include "gno.h"
#include "shapes.h"

/*!\brief Compute the pdf of the GNO distribution at a given location 
 * \param x position at which the pdf is computed
 * \param params parameter structure containing the parameters of the GNO
 * \return Value of the GNO at x, -1 in case of error (illegal shape parameter)
 */
double
pdfgno(double x, void *params)
{
    double y;
    sparam *p;
    
    p = (sparam *)params;

    if (p->alpha <= 0) {
      std::cout<<"Illegal value for distribution scale alpha: "<<
	p->alpha<<std::endl;
      return -1;
    }
    if (p->kappa == 0)
      y = (x - p->xi) / p->alpha;
    else
      y = -1. / p->kappa * log(1 - p->kappa * (x - p->xi) / p->alpha);
    
    return 1. / (p->alpha * sqrt(2 * M_PI)) * exp((p->kappa * y - y*y/2.));
}


/*!\brief Compute the cumulative distribution function of the GNO 
 * \param x position at which the cdf is computed
 * \param params parameter structure containing the parameters of the GNO
 * \return Value of the CDF at x, -1 in case of error (illegal shape parameter)
 */
double
cdfgno(double x, void *params)
{
    double a;   /* Lower integration limit if not -\infty */
    double cdf, error;
    gsl_function F;
    gsl_integration_workspace *w;
    sparam *p;
    
    p = (sparam *)params;

    if (p->alpha <= 0) {
      std::cout<<"Illegal value for distribution scale alpha: "<<
	p->alpha<<std::endl;
      return -1;
    }

    F.function = &pdfgno;
    F.params = params;
    w = gsl_integration_workspace_alloc(1024);
    if (p->kappa >= 0.0) {
	gsl_integration_qagil(&F, x, 1e-7, 1e-5, 1024, w, &cdf, &error);
	gsl_integration_workspace_free(w);
    } else {
	a = p->xi + p->alpha / p->kappa;
	gsl_integration_qags(&F, a, x, 1e-7, 1e-5, 1024, w, 
			    &cdf, &error);
    }
    gsl_integration_workspace_free(w);
    return cdf;
}


/*!\brief Compute the cumulative distribution function of the GNO minus a constant
 * \param x position at which the cdf is computed
 * \param params parameter structure containing the parameters of the GNO and the constant to be subtracted
 * \return Value of the CDF at x - params->zero, -99 in case of error (illegal shape parameter)
 */
double 
cdfgno_zero(double x, void *params)
{
    double val;
    sparam *p;

    p = (sparam *)params;
    val = cdfgno(x, params);
    if (val == -1)
	return -99;
    return val - p->zero;
}


/*!\brief Draw a random deviate from a generalized normal distribution
 * \param rng a GNU Scientific Library random number generator instance
 * \param parameters of the GNO
 * \return a random number, -99 in case of error, which sadly is not a unique indication of an error
 */
double
ran_gno(const gsl_rng *rng, void *params)
{
    /* This function draws random deviates from a generalized normal
       distribution. Since the GNO is neither analytically integrable
       nor log-concave, both the analytic transformation method and
       the adaptive rejection method cannot be used. Instead this
       function uses a numeric implementation of the transformation
       method. This is extremely slow and not really suited for
       massive Monte-Carlo problems. */
    int status;
    int iter, max_iter = 100;
    double r;
    double x_lo, x_hi;
    sparam *p;
    gsl_function F;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;

    p = (sparam *)params;
    F.function = &cdfgno_zero;
    F.params = params;
    if (p->kappa < 0) {
	x_lo = p->xi + p->alpha / p->kappa + 1e-9;
	x_hi = p->xi - 50 * p->kappa;
    } else {
	x_lo = p->xi - 50 * p->kappa;
	x_hi = p->xi + p->alpha / p->kappa - 1e-9;
    }

    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc(T);      
    /* Sometimes we request a random number so far out in the tail
       that it is outside our integration range. If this happens, draw
       another random deviate. */
    do {
	p->zero = gsl_rng_uniform(rng);
	status = gsl_root_fsolver_set(s, &F, x_lo, x_hi);
    } while (status == GSL_EINVAL);
    iter = 0;
    do {
	iter++;
	status = gsl_root_fsolver_iterate(s);
	if (status > 0) {
	  std::cout<<"Error in gsl_root_fsolver_iterate: "<<
	    gsl_strerror(status)<<std::endl;;
	  return -99;
	}
	r = gsl_root_fsolver_root(s);
	x_lo = gsl_root_fsolver_x_lower(s);
	x_hi = gsl_root_fsolver_x_upper(s);
	status = gsl_root_test_interval(x_lo, x_hi, 1e-7, 1e-3);
    } while(status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free(s);

    if (status == GSL_SUCCESS)
	return r;
    else
      std::cout<< "Could not obtain GNO random deviate"<<std::endl;
    return -99;
}
