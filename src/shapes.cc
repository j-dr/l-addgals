#include <ctime>
#include <istream>
#include <fstream>
#include <iostream>
#include "shapes.h"
#include "shapeconfig.h"
#include "calcparams.h"
#include "efunc.h"
#include "gno.h"

using namespace std;

istream & operator>>(istream & is, shapemag & in)
{
  is >> in.bands[0] >> in.bands[1] >> in.bands[2] >>
    in.bands[3] >> in.bands[4];
  return is;
}

//Given a vector of DES mags, assuming 6 mags per galaxy
//generate an ellipticity and shape for each.
void generate_shapes(vector<double> mag, vector<double> e, vector<double> s,
		     int nelem, int vl)
{
  int i, iam;
  gsl_rng **myrng;
  double mymag, e1, e2;
  prefstruct prefs;
  eparam eparams;
  sparam sparams;
  unsigned long int seed;


  //read in default prefs to a struct
  vector<string> vdp(dl);
  for (vector<string>::iterator itr=vdp.begin(); itr!=vdp.end(); itr++)
    {
      string str = string(default_prefs[(itr-vdp.begin())]);
      *itr = str;
    }

  parse_config(vdp, prefs);
  seed = static_cast<unsigned long int>(time(NULL));

  *myrng = (gsl_rng*)calloc(prefs.nthreads, sizeof(gsl_rng*));
  for(i = 0; i < prefs.nthreads; i++) {
      myrng[i] = gsl_rng_alloc(gsl_rng_ranlxd2);
      gsl_rng_set(myrng[i], seed + i);
  }

  //for now, nthreads always 1, fix iam
  iam = 0;

  //Might consider using OMP here?
  for (i=0; i<s.size(); i++)
    {
      mymag = mag[nelem - 1 + i * vl];
      calceparams(mymag, &eparams);
      calcsparams(mymag, &sparams);
      rng_efunc(myrng[iam], &eparams, &e1, &e2);
      e[2 * i] = e1;
      e[2 * i + 1] = e2;
      s[i] = ran_gno(myrng[iam], &sparams);
    }
  
  for(i = 0; i < prefs.nthreads; i++)
    gsl_rng_free(myrng[i]);

}

int main(int argc, char* argv[])
{
  int ngal, nelem, vl;

  nelem = 3;
  vl = 5;

  if (argc<2){
    cerr << "Usage: " << argv[0] << " magnitudes" << endl;
    return 1;
  }

  ifstream mag_file(argv[1]);
  vector<shapemag> smag;

  if (mag_file.fail()) {
    cerr<<"error: cannot open "<<argv[1]<<endl;
    exit(1);
  }

  copy(istream_iterator<shapemag>(mag_file),
       istream_iterator<shapemag>(),
       back_inserter(smag));

  ngal = smag.size();

  vector<double> e(2*ngal);
  vector<double> s(ngal);
  vector<double> mag(ngal*5);

  for (vector<double>::iterator itr=mag.begin(); itr!=mag.end(); itr++){
    *itr = smag[itr-mag.begin()].bands[(itr-mag.begin())%5];
  }

  generate_shapes(mag, e, s, nelem, vl);
  //write out magnitudes
  ofstream s_file("./des_stest.txt");
  ofstream e_file("./des_etest.txt");

  vector<double>::iterator ditr;

  for (ditr=s.begin(); ditr!=s.end(); ++ditr)
    {
      s_file << *ditr << "\n";
    }

  for (ditr=e.begin(); ditr!=e.end(); ++ditr)
    {
      e_file << *ditr;
      if ((ditr-e.begin()+1)%2==0)
	{
	  e_file << "\n";
	} else e_file << " ";
    }
}

