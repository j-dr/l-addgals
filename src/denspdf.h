//#include <iostream>
#include <string>
#include <vector>
#include <assert.h>
#include "interpolation.h"
//#include "fivetuple.h"
#include "stdafx.h"

using namespace alglib;

struct PDFentry{
  float mr;
  float params[5];
};

struct PDFarr{
  float z;
  std::vector <struct PDFentry> magbin;
};

class denspdf{
  public:

    void initalize(std::string file){
      //cout<<"Initializing the denspdf structure from file "<<file<<endl;
      read_file(file);
      make_splines();
    }

    float cm(float z, float mr){return spline2dcalc(spline_cm, z, mr);}
    float cs(float z, float mr){return spline2dcalc(spline_cs, z, mr);}
    float fm(float z, float mr){return spline2dcalc(spline_fm, z, mr);}
    float fs(float z, float mr){return spline2dcalc(spline_fs, z, mr);}
    float p(float z, float mr){return spline2dcalc(spline_p, z, mr);}

  private:
    std::vector <struct PDFarr> pdfhist;
    spline2dinterpolant spline_cm, spline_cs, spline_fm, spline_fs, spline_p;
    void make_splines();
    void read_file(std::string);

};
