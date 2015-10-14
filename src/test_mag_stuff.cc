#include <cmath>
#include <iostream>

static vector <double> magnitude;
static vector <double> lumnumdens;

void ReadLFFile(void){
  string filename = "../sdssinp/lumnumberdens.dat";
  ifstream file(filename.c_str());
  if (file.fail()) {
    cerr<<"error: cannot open "<<filename<<endl;
    exit(1);
  }
  else cout<<"reading "<<filename<<endl;
  int entries=3801;
  magnitude.resize(entries);
  lumnumdens.resize(entries);
  for (int ii=0; ii<entries; ii++) {
    file>>magnitude[ii]>>lumnumdens[ii];
  }
  file.close();
}

double LumNumberDensityInterp(double M){
  //if((M<-25)||(M>-5)){
  if((M<-25)||(M>1)){
    cerr<<"Magnitude out of bounds in LumNumberDensity function: "<<M<<endl;
    exit(1);
  }
  vector<double>::iterator pos;
  pos=upper_bound(magnitude.begin(), magnitude.end(), M);
  int i = distance(magnitude.begin(),pos);
  double delta_x = magnitude[i] - magnitude[i-1];
  double dx      = M        - magnitude[i-1];
  double delta_y = lumnumdens[i] - lumnumdens[i-1];
  double y = dx/delta_x*delta_y+lumnumdens[i-1];

  return y;
}


//Returns the magnitude of a specified number density
double NdensMagnitude(double ndens){
  double mmin = -23.;
  //double mmax = -5;
  double mmax = 0.;
  double diff = 1.;
  double nndensmid;
  double mmid;
  //cout<<ndens<<" "<<LumNumberDensity(mmin)<<" "<<LumNumberDensity(mmax)<<" ";
  while(diff>0.0001){
    mmid=(mmin+mmax)/2.;
    //    cout<<"range"<<ndens<<" "<<mmin<<" "<<mmax<<endl;
    double nndens1 = LumNumberDensityInterp(mmin);
    double nndens2 = LumNumberDensityInterp(mmax);
    nndensmid = LumNumberDensityInterp((mmin+mmax)/2);
    if((ndens<nndens1) || (ndens>nndens2)){
      cout<<mmin<<" "<<mmid<<" "<<mmax<<endl;
      cout<<"bad val"<<nndens1<<" "<<ndens<<" "<<nndens2<<endl;
    }
    assert((ndens>nndens1)&&(ndens<nndens2));
    diff = mmax-mmin;
    if(diff<0.0001) break;
    if(ndens>nndensmid) mmin=mmid;
    else mmax=mmid;
  }
  return mmid;
}

int main(void){
  ReadLFFiles();
  double mag = -19.;
  double phi = NdensMagnitude(phi);
  cout<<phi;
}
