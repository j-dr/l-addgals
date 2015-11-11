/*!\file calcparams.c
 * \brief Compute ellipticity pdf parameters
 */

#include <iostream>
#include "calcparams.h"
#include "efunc.h"
#include "gno.h"
#include "shapes.h"

/*!\brief Compute ellipticity pdf parameters
 * \param mag magnitude at which the pdf parameters are computed
 * \param p pointer to an ellipticity parameter structure
 * \return 0 if object in magnitude range, 1 if magnitude was reset to
 * minimum magnitude, 2 if magnitude was reset to maximum magnitude,
 * -1 if a computed parameter has an illegal value
 */
int
calceparams(double mag, eparam *p)
{
    int status = 0;
    float sigma_p;

    if (mag < prefs.eminmag) {
	mag = prefs.eminmag;
	status = 1;
    }
    if (mag > prefs.emaxmag) {
	mag = prefs.emaxmag;
	status = 2;
    }
    mag -= prefs.erefmag;

    p->b = prefs.p_0 + mag * (prefs.p_1 + mag * (prefs.p_2 + mag * prefs.p_3));
    sigma_p = prefs.sigma_0 + mag * (prefs.sigma_1 + mag * 
				     (prefs.sigma_2 + mag * prefs.sigma_3));
    if (sigma_p <= 0 || sigma_p > 3) {
      std::cout<<"Value for sigma_p" << sigma_p <<
	"outside allowed range at magnitude "<< mag << std::endl;
	status = -1;
    }
    p->a = sigma_p * pow(p->b, 1. / p->b);

    //printf("mag: %f sigma_p: %f p: %f\n", mag+prefs.erefmag, sigma_p, p->b);
    return status;
}


/*! \brief Compute size PDF parameters 
 * \param mag magnitude at which the pdf parameters are computed
 * \param p pointer to a size parameter structure
 * \return 0 if object in magnitude range, 1 if magnitude was reset to minimum magnitude, 2 if magnitude was reset to maximum magnitude, -1 if a computed parameter has an illegal value
 */
int
calcsparams(double mag, sparam *p)
{
    int status = 0;

    if (mag < prefs.sminmag) {
	mag = prefs.sminmag;
	status = 1;
    }
    if (mag > prefs.smaxmag) {
	mag = prefs.smaxmag;
	status = 2;
    }
    mag -= prefs.srefmag;

    p->xi = prefs.xi_0 + mag * (prefs.xi_1 + mag * (prefs.xi_2 + mag * 
						    (prefs.xi_3 + mag * 
						     prefs.xi_4)));
    p->alpha = prefs.alpha_0 + mag * (prefs.alpha_1 + mag * 
				      (prefs.alpha_2 + mag * 
				       (prefs.alpha_3 + mag * 
					prefs.alpha_4))); 
    if (p->alpha <= 0) {
      std::cout<<"Value for alpha "<<p->alpha<<
	" outside allowed range at magnitude "<<mag<<std::endl; 
	status = -1;
    }
    p->kappa = prefs.kappa_0 + mag * (prefs.kappa_1 + mag * 
				      (prefs.kappa_2 + mag * 
				       (prefs.kappa_3 + mag * 
					prefs.kappa_4))); 
    return status;
}
