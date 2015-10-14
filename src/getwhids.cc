#include <vector>
#include <iostream>
#include "particle.h"
#include "halo.h"
#include "biniostream.h"

void ReadHeader(binifstream& hpfile){

  float tmp=0;
  int tmpint=0;
  int record=0;
  hpfile>>record;
  for(int i=0;i<2;i++){
    hpfile>>tmp;
  }
  hpfile>>tmpint;
  hpfile>>tmp;
  hpfile>>tmp;
  hpfile>>record;
  hpfile>>record;
  int ng, ntot;
  hpfile>>ng>>ntot;
  hpfile>>record;
}  

void GetWHids(vector <Halo *> halos, vector <Particle *> particles){
  string hpinfoname;
  //if(box == W1a)
  //hpinfoname = Dir()+"hvd1.dat";
  //else if(box == W1b)
  //hpinfoname = Dir()+"hvd1e100.dat";
  //else if(box == W2)
  //hpinfoname = Dir()+"hvd2.dat";
  if(box == W384a)
    hpinfoname = Dir()+"eb14.r8.dat";
  else if(box == W384b)
    hpinfoname = Dir()+"eb15.r8.dat";
  else if(box == W384c)
    hpinfoname = Dir()+"eb16.r8.dat";
  else{
    cout<<"Inconsistent box/sim type"<<simulation<<"\t"<<simulation<<" "<<box<<endl;
    exit(0);
  }

  binifstream hpfile(hpinfoname.c_str());
  if (hpfile.fail()) {
    cerr<<"error: cannot open "<<hpinfoname<<endl;
    exit(0);
  }
  else cout<<"Reading "<<hpinfoname<<endl;
  string outfilename = Dir()+"ahid."+MakeString(box,1)+".1.dat";	
  ofstream pfile(outfilename.c_str());
  if (pfile.fail()) {
    cerr<<"error	: cannot open "<<outfilename<<endl;    exit(1);
  }
  else cout<<"Writing to "<<outfilename<<endl;
    //  hpfile.bigendian();
  ReadHeader(hpfile);
  int record,tmpint;
  hpfile>>record;
  cout<<"Reading for "<<record/4.<<" groups."<<endl;
  vector <int> np;
  np.reserve((int) (record/4.));
  np.push_back(0);
  for(int i=0;i<record/4.;i++){
    hpfile>>tmpint;  
    if (tmpint>nhpmin) {
      //cout<<tmpint<<endl;
      np.push_back(tmpint);
    }
  }
  hpfile>>record;
  cout<<record<<" "<<record/4.<<endl;
  cout<<"Read "<<np.size()<<" groups"<<endl;
  hpfile>>record;
  cout<<record<<" "<<record/4.<<endl;
  int pid=0;
  for(int hi=0;hi<halos.size();hi++){
    for(int j=0;j<np[hi];j++){
      hpfile>>pid;
      
      particles[pid]->Hid(hi);
      //if((j==0)||(j==1)){
      //cout<<pid<<" "
      //    <<particles[pid]->Distance(halos[hi]->Position())<<" "
      //    <<halos[hi]->R200()<<endl;
      //}
      //cout<<particles[pid]->Hid()<<" "<<particles[pid]->X()<<" "<<halos[hi]->X()<<" "
      //  <<particles[pid]->Y()<<" "<<halos[hi]->Y()<<" "<<endl;
	//	  <<halos[hi-1]->X()<<" "
	//<<halos[hi+1]->X()<<endl;
      //if(j==0){
	//particles[pid]->PosPrint();
      //      cout<<hi<<" "<<pid<<" "<<particles[pid]->Hid()<<endl;
      //}
    }
  }
  for(int pi=1;pi<particles.size();pi++){
    pfile<<particles[pi]->Hid()<<endl;
  }
}
