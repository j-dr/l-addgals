#ifndef _GNO_H_
#define _GNO_H_ 1

#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>

/*!
 * \struct sparam
 * \brief PDF parameters for size distribution
 */
struct sparam {
    double xi;     /*!< location */
    double alpha;  /*!< scale */
    double kappa;  /*!< shape */
    double zero;   /*!< zero for root finding */
};
typedef struct sparam sparam;

double pdfgno(double , void *);
double cdfgno(double , void *);
double cdfgno_zero(double , void *);
double ran_gno(const gsl_rng *, void *);

#endif
