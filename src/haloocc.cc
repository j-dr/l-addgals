#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>

#include "particle.h"
#include "halo.h"
#include "linreg.h"

static int runcount = 0;

using namespace std;

float mean(vector <int> x){
  float mean = 0.0;
  for(int i=0; i < x.size(); i++){
    mean += x[i];
  }
  mean = mean/x.size();
  return mean ;
}


float sqrtpairs(vector <int> x){
  float pairs = 0.0;
  for(int i=0; i < x.size(); i++){
    pairs += x[i]*(x[i]+1);
  }
  pairs = pairs/x.size();
  return sqrt(pairs);
}

//
void HaloOcc(vector <Halo *>& halos){
  //  string hoccfilename=out_path+"mngals.dat";
  //ofstream hoccfile(hoccfilename.c_str());
  int nbins = 10;
  float binsize=0.2;
  float binmin = 13.5;

  vector <float> binmid;
  vector <float> binmean;
  vector <float> binpairs;
  vector <float> binmeanms;
  vector <float> binmeanb;
  vector <float> binmeanm;
  vector <float> binmeand;
  int bin = nbins-1;
  int hi = 0;
  while(bin>=0){
    vector <int> ngals;
    vector <int> nms;
    vector <int> nbri;
    vector <int> nmid;
    vector <int> ndim;
    // cout<<"*** "<<bin<<" "<<hi<<" "<<ngals.size()<<endl;
    // while(hi<halos.size()){
    while(1){
      if(hi>=halos.size()){
	bin--;
	break;
      }
      Halo * h = halos[hi];
      hi++;
      int mybin =  (int) (floor((log10(h->M())-binmin)/binsize));
      //if(h->Ngal()){
      //cout<<"%% "<<bin<<" "<<mybin<<" "<<h->M()<<" "<<hi<<" "<<h->Ngal()<<" "
      //    <<h->Nbright()<<" "
      //    <<h->Ndim()<<" "<<endl;
      //}
      if(mybin==bin){
	if(h->Ngal()){
	  ngals.push_back(h->Ngal());
	  //	  hoccfile<<h->M()<<" "<<h->Ngal()<<endl;
	  //if(h->Nbright())
	  nms.push_back(h->Nmstar());
	  //if(h->Nmid())
	  nbri.push_back(h->Nbright());
	  //if(h->Nmid())
	  nmid.push_back(h->Nmid());
	  //if(h->Ndim())
	  ndim.push_back(h->Ndim());
	}
      }
      else if(mybin<bin){
	bin--;
	break;
      }
      else{
	//cout<<"Skipping high mass halo"<<log10(h->M())<<endl;
	break;
      }      
    }
    if(ngals.size()>0){
      //  cout<<binmid[bin+1]<<" "<<hi-1<<" "<<ngals.size()<<" "
      //	 <<mean(ngals)<<" "<<sqrtpairs(ngals)<<endl;
      binmid.push_back(binmin+(0.5+(bin+1))*binsize);
      binmean.push_back(mean(ngals));
      binmeanms.push_back(mean(nms));
      binmeanb.push_back(mean(nbri));
      binmeanm.push_back(mean(nmid));
      binmeand.push_back(mean(ndim));
      binpairs.push_back(sqrtpairs(ngals));
    }
  }
  LinearRegression lr;  
  LinearRegression lrb;  
  LinearRegression lrm;  
  LinearRegression lrd;  
  LinearRegression lrms;  
  //LinearRegression lrp;
  //  int runcount = 0;
#ifdef PARALLEL
  std::stringstream out_Task;
  out_Task<<ThisTask;
  string outstring = out_path+"hocc_fit."+MakeString(runcount, 3)+".dat"+out_Task.str();
#else
  string outstring = out_path+"hocc_fit."+MakeString(runcount, 3)+".dat";
#endif
  ofstream outfile(outstring.c_str());
#ifdef PARALLEL
  outstring = out_path+"hocc."+MakeString(runcount, 3)+".dat"+out_Task.str();
#else
  outstring = out_path+"hocc."+MakeString(runcount, 3)+".dat";
#endif
  ofstream outfile2(outstring.c_str());
  for(int i=0;i<binmid.size();i++){
    outfile2<<binmid[i]<<" "<<binmean[i]<<" "<<binpairs[i]<<" "<<binmeanms[i]<<" "<<binmeanb[i]<<" "<<binmeanm[i]<<" "<<binmeand[i]<<endl;
    lr.addPoint(Point2D(binmid[i],log10(binmean[i])));
    lrms.addPoint(Point2D(binmid[i],log10(binmeanms[i])));
    lrb.addPoint(Point2D(binmid[i],log10(binmeanb[i])));
    lrm.addPoint(Point2D(binmid[i],log10(binmeanm[i])));
    lrd.addPoint(Point2D(binmid[i],log10(binmeand[i])));
  }
  cout <<"### "<<pow(10.,-1*lr.getA()/lr.getB()) <<" "<<lr.getB()<<" "<<endl;
  outfile<<pow(10.,-1*lr.getA()/lr.getB()) <<" "<<lr.getB()<<" "<<endl;
  outfile<<pow(10.,-1*lrms.getA()/lrms.getB()) <<" "<<lrms.getB()<<" "<<endl;
  outfile<<pow(10.,-1*lrb.getA()/lrb.getB()) <<" "<<lrb.getB()<<" "<<endl;
  outfile<<pow(10.,-1*lrm.getA()/lrm.getB()) <<" "<<lrm.getB()<<" "<<endl;
  outfile<<pow(10.,-1*lrd.getA()/lrd.getB()) <<" "<<lrd.getB()<<" "<<endl;


  for(hi=0;hi<halos.size();hi++){

  }
  //  hoccfile.close();


}
  
