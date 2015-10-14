#include <iostream>
#include "biniostream.h"

using namespace std;

int main(void){
  int start = 0;
  int end = 0;
  const float LengthUnit = 3000.0;

  for(int i=start;i<=end;i++){
    for(int j=start;j<=end;j++){
      for(int k=start;k<=end;k++){
	char fstrh[20];
	char fstrhout[20];
	cout<<"[read_dens] looking for cube "<<i<<" "<<j<<" "<<k<<endl; 
	sprintf(fstrh,"/data2/evrard/HV/lcdm/POX/rnn/h.%02i.%02i.%02i",i,j,k);
	binifstream hfile(fstrh);
	short int tmp_dens;
	if (hfile.fail()) {
	  cerr<<"error: cannot open file '" <<fstrh<<"'"<<endl;
	  return 0;
	}
	else cout<<"reading file '"<<fstrh<<"'"<<endl;

	//sprintf(fstrhout,"/data2/risa/tmp_hvdens/h.%02i.%02i.%02i.txt",i,j,k);
	//ofstream hout(fstrhout);
	ofstream hout("ddd.out");
       	
	while(hfile){
	  //** Read 1 short int from density file
	  hfile>>tmp_dens;
	  float ihmax = 32767;
	  float ihmin = -ihmax;
	  float hmaxabs = 1.;
	  float hminabs = 1. / (16.*2.*ihmax);
	  float sfaci =  (log(hmaxabs) - log(hminabs))/  (ihmax - ihmin);
	  float dist = hminabs * exp( sfaci* (tmp_dens - ihmin) );
	  float dist8= dist*LengthUnit;
	
	  hout<<tmp_dens<<" "<<dist8<<endl;
	}
      }
    }
  }
  return 1;
}
	  
