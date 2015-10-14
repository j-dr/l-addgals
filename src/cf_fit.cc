#include "nr.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;
static  ofstream errout("errout.dat");

int cf_fit(string outfilename, ofstream &datfile, double &r0, double &gam, double &chi2){
  Vec_DP *xin_p, *yin_p, *sigin_p, *yn_p;
  //const int NPT=11; //fix this...
  const int NPT=7; //fix this...
  xin_p= new Vec_DP(NPT);
  yin_p= new Vec_DP(NPT);
  yn_p= new Vec_DP(NPT);
  sigin_p= new Vec_DP(NPT);
  Vec_DP &xin = *xin_p;
  Vec_DP &yin = *yin_p;
  Vec_DP &yn = *yn_p;
  Vec_DP &sigin = *sigin_p;


    //  Vec_DP xin(NPT), yin(NPT), sigin(NPT);
  int n0=0;
  ifstream infile(outfilename.c_str());
  //  cout<<"reading from "<<outfilename<<endl;
  if (infile.fail()) {
    cerr<<"[cf_fit]: cannot open "<<outfilename<<endl;    exit(1);
  }

  /*string tt = "cat "+outfilename;
  cout<<tt<<endl;
  system(tt.c_str());*/
  //LinearRegression lr;   // empty instance of linear regression
  //string command = "wc "+outfilename+" | awk '{print $1}' > wc.dat";
  //system(command);
  //ofstream wc("wc.dat");
  //wc>>NPT>>endl;
  
  for(int i=0;i<NPT;i++){
    infile>>xin[i]>>yin[i]>>sigin[i];
    // lr.addPoint(Point2D(log10(xin[i]),log10(yin[i])));
    //#ifdef DEBUG
    //cout<<xin[i]<<" "<<yin[i]<<" "<<sigin[i]<<endl;
    //#endif
    //  xin=r[i];
    //yin=xi[i];
    //sigin=err[i];
    if(yin[i]<=0) n0++;
  }
  
  //  cout <<"## "<<pow(10.,-1*lr.getA()/lr.getB()) <<" "<<lr.getB()<<" "<<lr.getStdErrorEst()<<endl;

  //  float central_slope = (log10(yin[4])-log10(yin[7]))/(log10(xin[4])-log10(xin[7]));
  float central_slope = (log10(yin[3])-log10(yin[6]))/(log10(xin[3])-log10(xin[6]));
  if(n0>0){
    	cout<<"negative values:"<<xin[4]<<" "<<yin[4]<<" "<<xin[7]<<" "<<yin[7]<<endl;
	//	datfile<<"-8 -8 -8 -8 -8 -8 -8";
	return 0;
  }
  

  /*
  int bad_type = 0;
  if(central_slope <-2.3){
    	errout<<"xi too steep:"<<xin[4]<<" "<<yin[4]<<" "<<xin[7]<<" "<<yin[7]<<endl;
	bad_type = -3;
	//return 0;
  }
  else if(central_slope >-1.6){
    	errout<<"xi too shallow:"<<xin[4]<<" "<<yin[4]<<" "<<xin[7]<<" "<<yin[7]<<endl;
	bad_type = -4;
	//return 0;
  }


  for(int i=1;i<NPT-1;i++){
    float g = (log10(yin[i-1])-log10(yin[i+1]))/(log10(xin[i-1])-log10(xin[i+1]));
    float xi_exp = pow(10.,(log10(yin[i-1])-g*(log10(xin[i-1])-log10(xin[i]))));
    float f = yin[i]/xi_exp;

    if(yin[i]>yin[i-1]){
      errout<<"non-decreasing cf:"<<xin[i-1]<<" "<<yin[i-1]<<" "<<xin[i]<<" "<<yin[i]<<endl;
      bad_type = -7;
    }
    if(i==7){
      if(yin[i]<1){
	errout<<"r0 too small:"<<xin[i-1]<<" "<<yin[i-1]<<" "<<xin[i]<<" "<<yin[i]<<endl;
	bad_type = -1;
      }
    }
    else if(i==10){
      if(yin[i]>1){
	errout<<"r0 too big:"<<xin[i-1]<<" "<<yin[i-1]<<" "<<xin[i]<<" "<<yin[i]<<endl;
	bad_type = -2;
      }
    }
    if(central_slope<-2.1){
      if((xin[i]>1)&&(xin[i]<10)){
        if((f>1.30)||(f<0.80)){
          errout<<"xi<-2.1 too bumpy:"<<xin[4]<<" "<<yin[4]<<" "<<xin[7]<<" "<<yin[7]<<endl;
	  bad_type = -5;
        }
      }
    }
    else{
      if(xin[i]<10){
        if((f>1.30)||(f<0.80)){
          errout<<"xi too bumpy:"<<central_slope<<" "<<xin[i]<<" "<<f<<" "<<g<<" "<<xi_exp<<" "<<yin[i]<<endl;
	  bad_type = -6;
        }
      }
    }
  }
  */

  //  const int NPT2 = NPT-n0;  
  /*  Vec_DP x(NPT2), y(NPT2), sig(NPT2);
  for(int i=0;i<NPT;i++){
    if(yin[i]>0){
      x[i] = log(xin[i]);
      sig[i] = sigin[i]/yin[i];
      y[i] = log(yin[i]);
    }
    else i--;
    }*/
  //  cout<<"fiting "<<NPT2<<" points"<<endl;
  bool mwt = true;
  double a=0, b=0, siga=0, sigb=0, q=0;
  //  NR::fit(x, y, sig, mwt, a, b, siga, sigb, chi2, q);
  for(int i = 0; i<NPT; i++){
      xin[i] = log(xin[i]);
      yn[i] = (log(yin[i]+sigin[i])+log(yin[i]-sigin[i]))/2.;
      sigin[i] = sigin[i]/yin[i];
      yin[i] = log(yin[i]);

  }
  NR::fit(xin, yin, sigin, mwt, a, b, siga, sigb, chi2, q);
  chi2 = chi2/(NPT-2);
  r0 = exp(-1*a/b);
  gam = b;

  float maxdev = 0;
  float maxdev2 = 0;
  for(int i = 0; i<NPT-1; i++){
    float dev1 = (exp(yin[i])/(pow((exp(xin[i])/r0),gam)));
    if(abs(dev1-1)>maxdev) maxdev = abs(dev1-1);
    float dev2 = (exp(yin[i])-(pow((exp(xin[i])/r0),gam)))/(sigin[i]*exp(yin[i]));
    if(abs(dev2-1)>maxdev2) maxdev2 = abs(dev2-1);
    //    cout<<exp(xin[i])<<" "<<exp(yin[i])<<" "<<pow((exp(xin[i])/r0),gam)<<" "<<dev<<" "<<abs(dev-1)<<endl;
  }
  //datfile<<setw(4);
  //datfile<<" "<<r0<<" "<<gam<<" "<<chi2<<" ";
    //<<central_slope<<" "<<exp(yin[1]-yin[6])/pow(exp(xin[1]-xin[6]), gam)<<" "
    // <<maxdev<<" "<<maxdev2<<" "<<q<<" "<<exp(yin[0])/(pow(exp(xin[0])/r0, gam))	 <<" "<<exp(yin[3])/(pow(exp(xin[3])/r0, gam))<<endl;

  //  cout<<exp(yin[3])/(pow(exp(xin[3])/r0, gam))<<" ";
  //  cout<<" "<<r0<<" "<<gam<<" "<<chi2<<" "<<exp(yin[1]-yin[6])/pow(exp(xin[1]-xin[6]), gam)<<" "<<maxdev<<" "<<maxdev2<<endl;
  //  NR::fit(xin, yn, sigin, mwt, a, b, siga, sigb, chi2, q);
  //r0 = exp(-1*a/b);
  //gam = b;
  //  datfile<<r0<<" "<<gam<<" "<<chi2/(NPT+2)<<" "<<q<<endl;

  //  cout<<"2 "<<r0<<" "<<gam<<" "<<chi2/(NPT2+2)<<endl;

  delete xin_p;
  delete yin_p;
  delete yn_p;
  delete sigin_p;


  return 1;
}
