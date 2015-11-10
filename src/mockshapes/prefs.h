#ifndef _PREFS_H_
#define _PREFS_H_ 1

#include <limits.h>

#include <fitsio.h>

typedef struct 
{
    char   magname[FLEN_KEYWORD];     /* Column name for unlensed magnitude */
    char   lensmagname[FLEN_KEYWORD]; /* lensed (magnified) magnitude name */
    int    nelem;           /* position of magnitude in vector */
    int    magonly;         /* set to true if catalog without shape and size */
    char   sname[FLEN_KEYWORD]; /* intrinsic size name in output */
    char   lenssname[FLEN_KEYWORD]; /* observed size name in output */
    double sminmag;         /* minimum magnitude of size poly fit */
    double srefmag;         /* reference magnitude of size poly fit */
    double smaxmag;         /* maximum magnitude of size poly fit */
    double inseeing;        /* seeing under which the model was determined */
    double outseeing;       /* seeing of the output catalog */
    double xi_0, xi_1, xi_2, xi_3, xi_4; /* polynomial parameters of xi */
    double alpha_0, alpha_1, alpha_2, alpha_3, alpha_4; /* polynomial parameters of alpha */
    double kappa_0, kappa_1, kappa_2, kappa_3, kappa_4; /* polynomial parameters of kappa */
    char   ename[FLEN_KEYWORD]; /* intrinsic ellipticity name in output */
    char   shearname[FLEN_KEYWORD]; /* intrinsic ellipticity name in output */
    double eminmag;         /* minimum magnitude of ellipticity poly fit */
    double erefmag;         /* reference magnitude of ellipticity poly fit */
    double emaxmag;         /* maximum magnitude of ellipticity poly fit */
    double epsmax;          /* maximum ellipticity */
    double p_0, p_1, p_2, p_3; /* polynomial parameters of b */
    double sigma_0, sigma_1, sigma_2, sigma_3; /* polynomial parameters of sigma */
    int   verbose;         /* level of verbosity: 0-5 */
    int   seed;            /* Random seed */
    int   nthreads;        /* Number of threads if compiled with OpenMP */
    char  incat[PATH_MAX]; /* Path to the input catalog */
    char  shearcat[PATH_MAX]; /* Path to the raytracing catalog */
} prefstruct;

prefstruct prefs;

/* Function in prefs.c */
void dumpprefs(void);
void readprefs(char *, char **, char **, int);
int findkeys(char *, char keyw[][16]);

#endif
