#include <vector>
#include "galaxy.h"
#include "halo.h"

void AddToHalos(vector <Galaxy *> &galaxies, vector <Halo *> &halos){
  for(int gi=0; gi<galaxies.size();gi++){
    int hi = galaxies[gi]->P()->Hid();
    //    cout<<hi<<endl;
    if(hi>=0){
      double halor = galaxies[gi]->P()->Distance(halos[hi]->Position());
      //if (galaxies[gi]->Central())
      //cout<<"a galaxy halor "<<halor<<endl;
      if(halor<halos[hi]->R200()){
#ifndef BCGS
	if(!halos[hi]->Setgals()){
	  halos[hi]->Brightest(gi);
	  halos[hi]->Central(gi);
	  halos[hi]->CentralDistance(halor);
	}
	else{
	  if(galaxies[gi]->Mr()<galaxies[halos[hi]->Brightest()]->Mr())
	    halos[hi]->Brightest(gi);
	  if(halor<halos[hi]->CentralDistance()){
	    halos[hi]->Central(gi);
	    halos[hi]->CentralDistance(halor);
	  }
	}
#endif
	//increment halo occupation statistics
	halos[hi]->IncrN(); 
	if(galaxies[gi]->Mr()<Mstar) halos[hi]->IncrNms(); 
	if(galaxies[gi]->Mr()<-20) halos[hi]->IncrNdim(); 
	if(galaxies[gi]->Mr()<-21) halos[hi]->IncrNmid(); 
	if(galaxies[gi]->Mr()<-22) halos[hi]->IncrNbri(); 
      }
    }
    else{
      //      cout<<"galaxy "<<gi<<" is not in a halo"<<endl;
    }
  }
}

void SwapGalaxies(vector <Galaxy *> &galaxies, vector <Halo *> &halos){
  cout<<"Adding to halos"<<endl;
  //AddToHalos(galaxies, halos);
  cout<<"Assigned to halos"<<endl;

  int halos_in_vol = 0;
  int halos_out_vol = 0;
  int central_galaxies =0;

  for(int hi=0; hi<halos.size();hi++){
    //    cout<<halos[hi]->Ra()<<" "<<halos[hi]->Dec()<<" "<<halos[hi]->Z();

    //	Galaxy * galaxy = new Galaxy(mag,ngal, ff, cm, cs);//, fms, fmi, fss, fsi);
    //galaxyslice.push_back(galaxies[gi]);

    if(halos[hi]->InVol()) halos_in_vol++;
    if(halos[hi]->Setgals()){
      if(!halos[hi]->InVol()){ 
	halos_out_vol++;
	//	cout<<halos[hi]->X()<<" "<<halos[hi]->Y()<<" "<<halos[hi]->Z()<<" "<<halos[hi]->Ra()<<" "<<halos[hi]->Dec()<<endl;
      }
#ifdef SWAP
      int cid = halos[hi]->Central();
      int bid = halos[hi]->Brightest();
      float cmr = galaxies[cid]->Mr();
      float bmr = galaxies[bid]->Mr();
      //cout<<halos[hi]->X()<<" "<<halos[hi]->Y()<<" "<<halos[hi]->Z()<<" "<<halos[hi]->Ra()<<" "<<halos[hi]->Dec()<<" g:";
      //cout<<cid<<" "<<galaxies[cid]->P()->X()<<" "<<galaxies[cid]->P()->Y()<<" "<<galaxies[cid]->P()->Z()<<" "
      //  <<galaxies[cid]->Ra()<<" "<<galaxies[cid]->Dec()<<" ";
      galaxies[cid]->Mr(bmr);
      galaxies[bid]->Mr(cmr);
      galaxies[cid]->Centralize(halos[hi]);
      //cout<<" after: "<<galaxies[cid]->P()->X()<<" "<<galaxies[cid]->P()->Y()<<" "<<galaxies[cid]->P()->Z()<<" "
      //  <<galaxies[cid]->Ra()<<" "<<galaxies[cid]->Dec()<<endl;
      central_galaxies++;
#endif
    }
#ifndef BCGS       // Kind of useless if using BCGs -- all appropriately massive halos should 
#ifndef SNAPSHOT // have a BCG, and if we use something like Millennium with halos down to  
                   // low mass (b/c of downsampling) there are lots of halos with no galaxies  
                   // so this just creates a lot of unnecessary chatter
    else if(halos[hi]->InVol()){
      cout<<hi<<" galaxies not set for halo "<<halos[hi]->M()<<endl;
    }
#endif
#endif
  }
  cout<<"halos:"<<halos_in_vol<<" "<<halos_out_vol<<" galaxies"<<central_galaxies<<" "<<endl;
}
