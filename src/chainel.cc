#include "chainel.h"

void ChainEl::Print()const{
    for(int i=0;i<NP;i++)
      cout<<setprecision(3)<<chainval[i]<<" ";
    cout<<"\t";
}

void ChainEl::Write(ofstream &outfile)const{
  for(int i=0;i<NP;i++)
    outfile<<setprecision(5)<<chainval[i]<<" ";
}

void ChainEl::PrintRotation()const{
  cout<<"rotation matrix:"<<endl;
  for(int i=0;i<NP;i++){
    for(int j=0;j<NP;j++){
      cout<<rotation_matrix[i][j]<<"\t";
    }
    cout<<endl;
  }
}

void ChainEl::SetStepSize(string covfile){
  cout<<"setting step size"<<endl;
  ifstream infile(covfile.c_str());
  if(!infile){
    cerr<<"cannot open "<<covfile<<endl;
    exit(1);
  }
  for(int i=0;i<NP;i++){
    for(int j=0;j<NP;j++){
      if(infile){
	infile>>rotation_matrix[i][j];
      }
      else {cout<<"matrix too small"<<endl;  exit(0);}
    }
  }
  cout<<"done setting step size"<<endl;
}

void ChainEl::Reset(ChainEl oldel){
  double step_mult[NP];
  double sig[NP];
  for(int i=0;i<NP;i++){
    // Gaussian step multiplied by the appropriate constant
    step_mult[i]=StepScaleFactor*NR::gasdev(cseed)*eigenvect[i];
  }
  for(int i=0;i<NP;i++){
    sig[i]=0.;
    for(int j=0;j<NP;j++){
      //	cout<<i<<" "<<j<<" "<<rotation_matrix[i][j]<<" "<<step_mult[j]<<" ";
      sig[i]+=rotation_matrix[i][j]*step_mult[j];
    }
    //      cout<<endl; 
  }
  for(int i=0;i<NP;i++){
    chainval[i] =  oldel.chainval[i] + sig[i]*use_param[i];
    //cout<<i<<" "<<chainval[i]<<" "<<oldel.chainval[i]<<" "<<sig[i]<<" "<<use_param[i]<<endl;
  }
}

int ChainEl::Okay()const{
    int okay = -1;
    for(int i=0;i<NP;i++){
      if(chainval[i] > param_max[i]) okay = i;
      if(chainval[i] < param_min[i]) okay = i;
    }
    return okay;
}
