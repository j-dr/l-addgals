#include "shapes.h"
#include <ctime>

using namespace std;

//Given a vector of DES mags, assuming 6 mags per galaxy
//generate an ellipticity and shape for each.
void generate_shapes(vector<double> mags, vector<double> e, vector<double> s)
{
  static char prefsname[4096];
  static prefstruct prefs = default_prefs;
  static int argc = 4;
  char** argv;

  int narg = 0;
  int i, a, iam, vl;
  gsl_rng **myrng;
  double mymag, e1, e2;
  eparam eparams;
  sparam sparams;
  
  //kludge to mess with as little of mockshapes as necessary
  for (i=0;i<argc;i++){
    argv[i] = calloc(4096, sizeof(char));
  }

  strcpy(argv[0], "mockshapes");
  strcpy(argv[1], "-SEED");
  strcpy(argv[2], asctime(localtime(&time(nullptr))));
  strcpy(argv[3], "pass");
    
  strcpy(prefsname, "mockshapes_default.conf");

  narg = 0;
  for(a = 1; a < argc; a++) {
      if (*(argv[a]) == '-') {
          opt = (int) argv[a][1];
          if (strlen(argv[a]) < 3 || opt == '-') {
              if (opt == '-')
                  opt = (int) tolower((int) argv[a][2]);
              switch (opt) {
              case 'c':
                  if (a < (argc - 1))
                      strcpy(prefsname, argv[++a]);
                  break;
              case 'd':
                  dumpprefs();
                  exit(EXIT_SUCCESS);
                  break;
              case 'h':
              default:
                  msg_message(MSG_FATAL, 0, "SYNTAX:\n%s", SYNTAX);
              }
          } else {
              argkey[narg] = &argv[a][1];
              argval[narg++] = argv[++a];
          }
      } else {
          snprintf(prefs.incat, PATH_MAX, "%s", argv[a]);
      }
  }

  //read in default prefs to a struct
  readprefs(prefsname, argkey, argval, narg);

  myrng = calloc(prefs.nthreads, sizeof(gsl_rng));
  for(i = 0; i < prefs.nthreads; i++) {
      myrng[i] = gsl_rng_alloc(gsl_rng_ranlxd2);
      gsl_rng_set(myrng[i], prefs.seed + i);
  }

  //for now, nthreads always 1, fix iam
  iam = 0;
  //using des mags, 6 mags per row
  vl = 6;

  //Might consider using OMP here?
  for (i=0; i<s.size(); i++)
    {
      mymag = mag[prefs.nelem - 1 + i * vl];
      calceparams(mymag, &eparams);
      calcsparams(mymag, &sparams);
      rng_efunc(myrng[iam], &eparams, &e1, &e2);
      e[2 * i] = e1;
      e[2 * i + 1] = e2;
      s[i] = ran_gno(myrng[iam], &sparams);
    }
}

