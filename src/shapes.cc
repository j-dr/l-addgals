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
void generate_shapes(vector<float>& mag, vector<double>& e, vector<double>& s,
		     int nelem, int vl)
{
  cout << "generating shapes" << endl;  
  int i, iam;
  gsl_rng *myrng;
  double mymag, e1, e2;
  prefstruct prefs;
  eparam eparams;
  sparam sparams;

  cout << "reading in default prefs" << endl;
  //read in default prefs to a struct
  vector<string> vdp(dl);
  for (vector<string>::iterator itr=vdp.begin(); itr!=vdp.end(); itr++)
    {
      string str = string(default_prefs[(itr-vdp.begin())]);
      *itr = str;
    }
  
  parse_config(vdp, prefs);

  long int seed  = (long int)time(NULL);
  
  //myrng = (gsl_rng**)malloc(prefs.nthreads*sizeof(gsl_rng*));
  myrng = gsl_rng_alloc(gsl_rng_ranlxd2);
  gsl_rng_set(myrng, seed);

#ifdef DEBUG_SHAPES
  ofstream s_file("./stest.txt");
  ofstream e_file("./etest.txt");
#endif

  if (!myrng) {
    cout << "Null pointer!"<<endl;
    exit(1);
  }
  /*cout << "instantiating RNG" << endl;
  for(i = 0; i < prefs.nthreads; i++) {
      myrng[i] = gsl_rng_alloc(gsl_rng_ranlxd2);
      gsl_rng_set(myrng[i], seed + i);
      }*/

  //for now, nthreads always 1, fix iam
  iam = 0;


  //Might consider using OMP here?
  for (i=0; i<s.size(); i++)
    {
      mymag = mag[nelem - 1 + i * vl];
      if (i==0) cout << "calculating e params" << endl;
      if (isnan(mymag))
	{
	  e[2 * i] = -99;
	  e[2 * i + 1] = -99;
	  continue;
	}
      calceparams(mymag, &eparams, prefs);
#ifdef DEBUG_SHAPES
      e_file << eparams.a << " " << eparams.b << " ";
      e_file.flush();
#endif
      if (i==0) cout << "drawing e" << endl;
      rng_efunc(myrng, &eparams, &e1, &e2, prefs);
#ifdef DEBUG_SHAPES
      e_file << e1 << " " << e2 << endl;
#endif

      e[2 * i] = e1;
      e[2 * i + 1] = e2;
    }

  for (i=0; i<s.size(); i++)
    {
      mymag = mag[nelem - 1 + i * vl];
      if (i==0) cout << "calculating s params" << endl;
      if (isnan(mymag))
	{
	  s[i] = -99;
	  continue;
	}
      calcsparams(mymag, &sparams, prefs);
#ifdef DEBUG_SHAPES
      s_file << sparams.xi << " " << sparams.alpha << " "
	     << sparams.kappa << " " << sparams.zero << " ";
      s_file.flush();
#endif
      if (i==0) cout << "drawing s" << endl;
      s[i] = ran_gno(myrng, &sparams);
#ifdef DEBUG_SHAPES
      s_file << s[i] << endl;
#endif      
    }
  
  cout << "freeing RNG" << endl;
  /*for(i = 0; i < prefs.nthreads; i++)
    gsl_rng_free(myrng[i]);*/
  gsl_rng_free(myrng);

}

#ifndef BCC
int main(int argc, char* argv[])
{
  int ngal, nelem, vl, i;

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

  cout<<"reading in mags"<<endl;

  copy(istream_iterator<shapemag>(mag_file),
       istream_iterator<shapemag>(),
       back_inserter(smag));

  ngal = smag.size();

  vector<double> e(2*ngal);
  vector<double> s(ngal);
  vector<float> mag(ngal*5);

  cout<< "copying shapemag struct into double array" << endl;
  for (vector<float>::iterator itr=mag.begin(); itr!=mag.end(); itr++){
    *itr = smag[(itr-mag.begin())/5].bands[(itr-mag.begin())%5];
  }

  cout << "First few mags: ";
  for (i=0; i<10; i++){
    cout<<mag[i]<<" ";
    if ((i+1)%5==0) cout << endl;
  }

  generate_shapes(mag, e, s, nelem, vl);
  
  //write out shapes
  cout << "writing shapes" << endl;
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
#endif
