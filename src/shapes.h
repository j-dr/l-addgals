#ifndef SHAPES_H
#define SHAPES_H

#include <vector>

#define FLEN_KEYWORD 2048

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

char *default_prefs[] = {
    "# Default configuration file for mockshapes",
    "# JPD 12.01.2010",
    "#",
    "#---------------------------- Input Catalog ------------------------------------",
    "MAGNITUDE     TMAG        # Column name for magnitude",
    "LENSMAGNITUDE LMAG        # Column name for lensed magnitude",
    "NMAG          3           # Position of magnitude in vector",
    "MAGNITUDE_ONLY N          # Set to Y if no shapes/sizes in catalog",
    " ",
    "#---------------------------- Size Parameters ----------------------------------",
    "SIZE          TSIZE       # Name of intrinsic (true) size column in output", 
    "LENSSIZE      OSIZE       # Name of observed (magnified) object size",
    "SMINMAG       17.5        # Minimum valid magnitude of polynomial fit",
    "SREFMAG       22.0        # Reference magnitude of polynomial fit",
    "SMAXMAG       26.24       # Maximum valid magnitude of polynomial fit",
    "INSEEING      0.635       # Seeing in which the model parameters were measured",
    "OUTSEEING     0.9         # Seeing of the output catalog",
    "XI_0          0.676538    # Parameters for 4th order polynomial fit to xi",
    "XI_1          -0.0597809  #",
    "XI_2          -0.00242377 #",
    "XI_3          -0.00276862 #",
    "XI_4          0.000889298 #",
    "ALPHA_0       0.175922    # Parameters for 4th order polynomial fit to alpha",
    "ALPHA_1       -0.0172894  #",
    "ALPHA_2       -0.00162605 #",
    "ALPHA_3       -0.00233168 #",
    "ALPHA_4       0.000617415 #",
    "KAPPA_0       -0.658123   # Parameters for 4th order polynomial fit to kapp",
    "KAPPA_1       -0.00185666 #",
    "KAPPA_2       0.0318665   #",
    "KAPPA_3       0.0010374   #",
    "KAPPA_4       -0.00132519 #",
    " ",
    "#--------------------------- Shape Parameters ----------------------------------",
    "ELLIPTICITY   TE          # Name of ellipticity column in output",
    "SHEAR         EPSILON     # Name of the observed shear estimator",
    "EMINMAG       21.38       # Minimum valid magnitude of polynomial fit",
    "EREFMAG       24.0        # Reference magnitude of polynomial fit",
    "EMAXMAG       26.81       # Maximum valid magnitude of polynomial fit",
    "EPSMAX        1.0         # Maximum ellipticity",
    "SIGMA_0       0.2679      # Parameters for cubic polynomial fit to sigma_p",
    "SIGMA_1       0.0415      #",
    "SIGMA_2       0.0010      #",
    "SIGMA_3       -0.0010     #",
    "P_0           1.2565      # Parameters for cubic polynomial fit to p",
    "P_1           0.0937      #",
    "P_2           -0.0049     #",
    "P_3           -0.0029     #",
    " ",
    "#---------------------------- Program Setup ------------------------------------",
    "VERBOSE       4           # Verbosity of output",
    "SEED          0           # Random number seed",
    "NTHREADS      0           # Number of OpenMP threads (0: use system default)",
    ""
};

prefstruct prefs;

void generate_shapes(std::vector<double> mags, std::vector<double> e, std::vector<double> s);

#endif
