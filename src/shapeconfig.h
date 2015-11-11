#ifndef SHAPECONFIG_H
#define SHAPECONFIG_H

#include <vector>
#include <string>

typedef struct 
{
    int    nelem;           /* position of magnitude in vector */
    double sminmag;         /* minimum magnitude of size poly fit */
    double srefmag;         /* reference magnitude of size poly fit */
    double smaxmag;         /* maximum magnitude of size poly fit */
    double inseeing;        /* seeing under which the model was determined */
    double outseeing;       /* seeing of the output catalog */
    double xi_0, xi_1, xi_2, xi_3, xi_4; /* polynomial parameters of xi */
    double alpha_0, alpha_1, alpha_2, alpha_3, alpha_4; /* polynomial parameters of alpha */
    double kappa_0, kappa_1, kappa_2, kappa_3, kappa_4; /* polynomial parameters of kappa */
    double eminmag;         /* minimum magnitude of ellipticity poly fit */
    double erefmag;         /* reference magnitude of ellipticity poly fit */
    double emaxmag;         /* maximum magnitude of ellipticity poly fit */
    double epsmax;          /* maximum ellipticity */
    double p_0, p_1, p_2, p_3; /* polynomial parameters of b */
    double sigma_0, sigma_1, sigma_2, sigma_3; /* polynomial parameters of sigma */
    int   seed;            /* Random seed */
    int   nthreads;        /* Number of threads if compiled with OpenMP */
} prefstruct;

extern int dl;
extern char* default_prefs[];

void parse_config(std::vector<std::string> config, prefstruct& out);

#endif
