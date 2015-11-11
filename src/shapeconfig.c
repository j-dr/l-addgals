#include <iostream>
#include <iomanip>
#include <string>
#include <exception>
#include <fstream>
#include <vector>
#include <algorithm>
#include "shapeconfig.h"

int dl = 54;
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

void parse_config(std::vector<std::string> config, prefstruct& out) {
  std::string str;

  for (std::vector<std::string>::iterator itr=config.begin(); itr!=config.end(); itr++)
    {
      str = *itr;
      std::string::size_type begin = str.find_first_not_of(" \f\t\v");
      //Skips blank lines
      if(begin == std::string::npos)
	continue;
      //Skips #
      if(std::string("#").find(str[begin]) != std::string::npos)
	continue;

      std::string firstWord;
      std::string secondWord;
      size_t kend = str.find(" ");
      try {
	firstWord = str.substr(0,str.find(" "));
      }
      catch(std::exception& e) {
	firstWord = str.erase(str.find_first_of(" "),str.find_first_not_of(" "));
      }
      secondWord = str.substr(kend+1);

      //remove whitespace
      firstWord.erase(std::remove_if(firstWord.begin(), firstWord.end(), isspace), firstWord.end());
      secondWord.erase(std::remove_if(secondWord.begin(), secondWord.end(), isspace), secondWord.end());
      std::transform(firstWord.begin(),firstWord.end(),firstWord.begin(), ::toupper);

      if(firstWord == "NELEM")
	out.nelem = atoi(secondWord.c_str());
      if(firstWord == "SMINMAG")
	out.sminmag = atof(secondWord.c_str());
      if(firstWord == "SREFMAG")
	out.srefmag = atof(secondWord.c_str());
      if(firstWord == "SMAXMAG")
	out.smaxmag = atof(secondWord.c_str());
      if(firstWord == "INSEEING")
	out.inseeing = atof(secondWord.c_str());
      if(firstWord == "OUTSEEING")
	out.outseeing = atof(secondWord.c_str());
      if(firstWord == "XI_0")
	out.xi_0 = atof(secondWord.c_str());
      if(firstWord == "XI_1")
	out.xi_1 = atof(secondWord.c_str());
      if(firstWord == "XI_2")
	out.xi_2 = atof(secondWord.c_str());
      if(firstWord == "XI_3")
	out.xi_3 = atof(secondWord.c_str());
      if(firstWord == "XI_4")
	out.xi_4 = atof(secondWord.c_str());
      if(firstWord == "ALPHA_0")
	out.alpha_0 = atof(secondWord.c_str());
      if(firstWord == "ALPHA_1")
	out.alpha_1 = atof(secondWord.c_str());
      if(firstWord == "ALPHA_2")
	out.alpha_2 = atof(secondWord.c_str());
      if(firstWord == "ALPHA_3")
	out.alpha_3 = atof(secondWord.c_str());
      if(firstWord == "ALPHA_4")
	out.alpha_4 = atof(secondWord.c_str());
      if(firstWord == "KAPPA_0")
	out.kappa_0 = atof(secondWord.c_str());
      if(firstWord == "KAPPA_1")
	out.kappa_1 = atof(secondWord.c_str());
      if(firstWord == "KAPPA_2")
	out.kappa_2 = atof(secondWord.c_str());
      if(firstWord == "KAPPA_3")
	out.kappa_3 = atof(secondWord.c_str());
      if(firstWord == "KAPPA_4")
	out.kappa_4 = atof(secondWord.c_str());
      if(firstWord == "EMINMAG")
	out.eminmag = atof(secondWord.c_str());
      if(firstWord == "EREFMAG")
	out.erefmag = atof(secondWord.c_str());
      if(firstWord == "EMAXMAG")
	out.emaxmag = atof(secondWord.c_str());
      if(firstWord == "EPSMAX")
	out.epsmax = atof(secondWord.c_str());
      if(firstWord == "P_0")
	out.p_0 = atof(secondWord.c_str());
      if(firstWord == "P_1")
	out.p_1 = atof(secondWord.c_str());
      if(firstWord == "P_2")
	out.p_2 = atof(secondWord.c_str());
      if(firstWord == "P_3")
	out.p_3 = atof(secondWord.c_str());
      if(firstWord == "SIGMA_0")
	out.sigma_0 = atof(secondWord.c_str());
      if(firstWord == "SIGMA_1")
	out.sigma_1 = atof(secondWord.c_str());
      if(firstWord == "SIGMA_2")
	out.sigma_2 = atof(secondWord.c_str());
      if(firstWord == "SIGMA_3")
	out.sigma_3 = atof(secondWord.c_str());
      if(firstWord == "SEED")
	out.seed = atof(secondWord.c_str());
      if(firstWord == "NTHREADS")
	out.nthreads = atof(secondWord.c_str());
    }
}
