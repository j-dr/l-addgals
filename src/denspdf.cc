#include <iostream>
#include <fstream>
#include <vector>
#include "denspdf.h"
#include "interpolation.h"

using namespace std;

//sets up the spline for general use
void denspdf::make_splines(void){

  real_1d_array zarr, mrarr;
  real_1d_array pdf_cm, pdf_cs, pdf_fm, pdf_fs, pdf_p;

  //fill our arrays with the appropriate coordinate information
  int nz = pdfhist.size();
  int nmag = pdfhist[0].magbin.size();
  cout<<nz<<" "<<nmag<<endl;
  zarr.setlength(nz);
  mrarr.setlength(nmag);
  pdf_cm.setlength(nz*nmag);
  pdf_cs.setlength(nz*nmag);
  pdf_fm.setlength(nz*nmag);
  pdf_fs.setlength(nz*nmag);
  pdf_p.setlength(nz*nmag);
  for (int i=0;i<nz;i++) zarr[i] = pdfhist[i].z;
  for (int i=0;i<nmag;i++) mrarr[i] = pdfhist[0].magbin[i].mr;

  //save the PDF values in arrays...
  int ind = 0;
  for (int j=0;j<nmag;j++){
    for(int i=0; i<nz; i++){
      pdf_cm[ind] = pdfhist[i].magbin[j].params[0];
      pdf_cs[ind] = pdfhist[i].magbin[j].params[1];
      pdf_fm[ind] = pdfhist[i].magbin[j].params[2];
      pdf_fs[ind] = pdfhist[i].magbin[j].params[3];
      pdf_p[ind] = pdfhist[i].magbin[j].params[4];
      //if (mrarr[j] == -20.0) cout<<pdfhist[i].magbin[j].params[0]<<endl;
      ind++;
    }
  }  

  //build our cubic splines
  spline2dbuildbicubicv(zarr, nz, mrarr, nmag, pdf_cm, 1, spline_cm);
  spline2dbuildbicubicv(zarr, nz, mrarr, nmag, pdf_cs, 1, spline_cs); 
  spline2dbuildbicubicv(zarr, nz, mrarr, nmag, pdf_fm, 1, spline_fm); 
  spline2dbuildbicubicv(zarr, nz, mrarr, nmag, pdf_fs, 1, spline_fs); 
  spline2dbuildbicubicv(zarr, nz, mrarr, nmag, pdf_p, 1, spline_p); 

}

//reads in the denspdf file, creates the histogram array
void denspdf::read_file(std::string file){

  int nz, nmag;
  float this_z, this_mag;
  float tparams[5];
  int i,j;

  //open the file to read
  const char *fcharp = file.c_str();
  ifstream inFile(fcharp, ios::in | ios::binary);  
  if(!(inFile.is_open())){
    cout<<"Error!  Couldn't find denspdf file "<<file<<endl;
    cout<<"Aborting."<<endl;
    exit;
  }

  //read the number of redshifts, allocate the pdf array
  inFile.read((char *)&nz,sizeof(int));
  pdfhist.resize(nz);

  //read in the redshifts, save to the arrays
  for(i=0;i<nz;i++){
    inFile.read((char *)&this_z,sizeof(float));
    pdfhist[i].z = this_z;
  }

  //read in the magnitude array information
  inFile.read((char*)&nmag, sizeof(int));
  for(i=0;i<nz;i++){
    pdfhist[i].magbin.resize(nmag);
    for(j=0;j<nmag;j++){
      if (i == 0){
        inFile.read((char*)&this_mag,sizeof(float));
      } else {
        this_mag = pdfhist[0].magbin[j].mr;
      }
      pdfhist[i].magbin[j].mr = this_mag;
    }
  }

  //finally we read in the pdf parameters
  for(i=0;i<nz;i++){
    for(j=0;j<nmag;j++){
      inFile.read((char*)&tparams,5*sizeof(float));
      for(int k=0;k<5;k++) pdfhist[i].magbin[j].params[k] = tparams[k];
    }
  }

  inFile.close();
}

