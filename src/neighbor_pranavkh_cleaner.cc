#include <iostream>
#include <fstream>
#include <vector>
//#include <omp.h>
#include <ANN/ANN.h>			// ANN declarations
#include "nr.h"
#include "galaxy.h"
#include "singleton.h"
//#include "choose.h"
#include <math.h>

int findCloseGalaxies2(vector <GalSED> &v, float mag, float den, float ThisZ, int ThisBCGs);


vector <GalSED> ReadSED();

#ifdef COLORS_FROM_RELATIVE_DENSITY
class GalaxyPercent{
public: 
  GalaxyPercent(float nndist, int gid):nndist(nndist), id(gid){};
  float dens()const{return nndist;};
  int gid()const{return id;};
  void dens(float tnndist){
    nndist = tnndist;
  }
  void gid(int tid){
    id = tid;
  }
private:
  float nndist;
  int id;
};

int CompareByDens (GalaxyPercent a, GalaxyPercent b)
{
  //return a->dens() > b->dens();
  return a.dens() > b.dens();
}

vector <float> GetNeighborPercents(vector <float> nndist, vector <Galaxy *> &galaxies)
{
  vector <float> color_percent(galaxies.size());
  int nColorBins = floor(zmax->GetVal()/ColorBinSize)+1;
  struct GalaxyPercent tGalaxyId(0.,0);
  vector <GalaxyPercent> GalaxyId;

  int nInColorBins = 0;
  int nInColorBin[nColorBins];
  for (int i=0;i<nColorBins;i++)
    nInColorBin[i] = 0;
  int min_gid = 1000000;
  int max_gid = 0;
  int min_i = 1000000;
  int max_i = 0;

  cout<<" Looping through percent bins."<<endl;
  for (int bin=0;bin<nColorBins;bin++){
    float binmin = bin*ColorBinSize;
    float binmax = binmin + ColorBinSize;

    for (int i=0;i<galaxies.size();i++){
      if (galaxies[i]->Z() < binmin || galaxies[i]->Z() > binmax)
	continue;
      float tnndist = nndist[i];
      tGalaxyId.dens(tnndist);
      tGalaxyId.gid(i);
      GalaxyId.push_back(tGalaxyId);
      nInColorBin[bin]++;
    }

    //check how many 0's are in the dens list
    int nZero = 0;
    for (int i=0;i<nInColorBin[bin];i++)
      if (GalaxyId[i].dens() <= 0.)
	nZero++;
    if (nZero > 0)
      cout<<" Before sorting there are "<<nZero<<" 0's in the anticipated range."<<endl;

    sort(GalaxyId.begin(),GalaxyId.end(),CompareByDens);

    //check how many 0's are in the dens list
    nZero = 0;
    for (int i=0;i<nInColorBin[bin];i++)
      if (GalaxyId[i].dens() <= 0.)
	nZero++;
    if (nZero > 0)
      cout<<" After sorting there are "<<nZero<<" 0's in the anticipated range."<<endl;


    for (int i=0;i<nInColorBin[bin];i++){
      color_percent[GalaxyId[i].gid()] = float(i)/float(nInColorBin[bin]);
      if (GalaxyId[i].gid() > max_gid)
	max_gid = GalaxyId[i].gid();
      if (GalaxyId[i].gid() < min_gid)
	min_gid = GalaxyId[i].gid();
    }

    GalaxyId.clear();

    cout<<"     This bin had "<<nInColorBin[bin]<<" particles in it."<<endl;
    nInColorBins += nInColorBin[bin];
    //cout<<"  done with the bin."<<endl;
  }
  cout<<"First 3 color percents: "<<color_percent[0]<<" "<<color_percent[1]<<" "<<color_percent[2]<<endl;
  cout<<"Last 3 color percents: "<<color_percent[galaxies.size()-1]<<" "<<color_percent[galaxies.size()-2]<<" "<<color_percent[galaxies.size()-1]<<endl;
  cout<<"Did the percent for a total of "<<nInColorBins<<" Galaxies (tot # = "<<galaxies.size()<<")"<<endl;
  cout<<"Range of galaxy id's covered: "<<min_gid<<" - "<<max_gid<<endl;
  return color_percent;
}
#endif

// new 2d version
// this calculates the proper ra-dec distance
// Galaxies:  The n nearest objects to the points
// Points:  The locations from which the distances are measured


vector <float> GetNeighborDist(vector <Galaxy *> galaxies, vector <Galaxy *> points){
  cout<<"In GetNeighborDist"<<endl;
  cout<<"  number of galaxies: "<<galaxies.size()<<endl;;
  cout<<"  number of points:   "<<points.size()<<endl;
  //PRNTV(galaxies.size());
#ifdef SNAPSHOT
  //in case of snapshot we project along the z-axis
  int nbins = (int) (ceil(sim.Boxsize()/zbsize));  
#else
  PRNTV(zmax->GetVal());
  int nbins = (int) (ceil(zmax->GetVal()/zbsize));
#endif
  PRNTV(nbins);
  vector <int> ginbin(nbins+1);
  for(int i=0;i<ginbin.size();i++){
    ginbin[i] = 0;
  }
  int OutNext = 0;
  cout<<" initial ginbin's: "<<ginbin[0]<<" "<<ginbin[1]<<" "<<ginbin[2]<<" "<<ginbin[3]<<" "<<ginbin[4]<<endl;
  for(int i=0;i<galaxies.size();i++){
    if(i>OutNext){
      //cout<<i<<" of "<<galaxies.size()<<endl;;
      OutNext += 1000;
    }
    ginbin[galaxies[i]->Zbin()]++;
    if (i < 5){
      cout<<galaxies[i]->Ra()<<" "<<galaxies[i]->Dec()<<" "<<galaxies[i]->Z()<<" "<<galaxies[i]->Zbin()<<" "<<ginbin[galaxies[i]->Zbin()]<<endl;
    }
  }
  int ginanybin = 0;
  for(int i=0;i<ginbin.size();i++){
    ginanybin += ginbin[i];
    //cout<<"  "<<i<<" "<<ginbin[i]<<endl;
  }
  cout<<" Number of galaxies in bins (should only be ngals if Magmin >= Magmin_dens): "<<ginanybin<<endl;

  cout<<"creating an int vector of size "<<nbins+1<<endl;
  vector <int> pinbin(nbins+1);
  cout<<"Zeroing the vector."<<endl;
  for(int i=0;i<pinbin.size();i++){
    pinbin[i] = 0;
  }
  cout<<"  Zeroed pinbin."<<endl;
  OutNext = 0;
  for(int i=0;i<points.size();i++){
    //cout<<i<<" of "<<galaxies.size()<<endl;
    if(i>OutNext){
      //cout<<i<<" of "<<galaxies.size()<<endl;;
      OutNext += 1000;
    }
    pinbin[points[i]->Zbin()]++;
    if (i < 5){
      cout<<points[i]->Ra()<<" "<<points[i]->Dec()<<" "<<points[i]->Z()<<" "<<points[i]->Zbin()<<" "<<pinbin[points[i]->Zbin()]<<endl;
    }
  }
  cout<<"  Set pinbin"<<endl;
  int pinanybin = 0;
  for(int i=0;i<pinbin.size();i++)
    pinanybin += pinbin[i];
  cout<<" Number of points in bins: "<<pinanybin<<endl;

#ifdef DEBUG
  cout<<"[neighbor] "<<nbins<<" z bins"<<endl;
#endif
  int startbin=0;
  bool have_started = false;
  for(int i=0;i<nbins;i++){
#ifdef DEBUG
        cout<<"[neighbor] "<<i<<" "<<ginbin[i]<<endl;
#endif
	//	if((ginbin[i]<10)&&(i>1)){
	if(ginbin[i]<10){
	  if(have_started == false){
	    startbin++;
	    //cout<<"not enough objects in bin "<<i<<endl;
	  }
	  else{
	    //cout<<"not enough objects in bin "<<i<<endl;
	    nbins = i;
	    break;
	  }
	}
	else have_started = true;
  }
  cout<<"using "<<startbin<<" "<<nbins<<endl;
  //vector <float> nndist(galaxies.size());
  vector <float> nndist(points.size());
  ANNpoint            query_pt;               // query point
  ANNidxArray         nn_idx;                 // near neighbor indices
  ANNdistArray        dists;                  // near neighbor distances
  ANNkd_tree * the_tree;
  double eps=0;
  int dim = 2;
  int k_want = 10;  // 10th nearest neighbor (first one is self)
  if(dens_measure == FIFTH)
    k_want = 5;

  query_pt = annAllocPt(dim);             // allocate query point
  //  int jmax = 0;
  string fnddout = out_path+"dd2.out";
  string fnnfout = out_path+"nf.out";
  ofstream ddout(fnddout.c_str());
  ofstream nfout(fnnfout.c_str());
  //loop over the z slices
  int  pi=0;
#ifdef DEBUG
  cout<<"Looping over "<<nbins-startbin<<"z bins"<<endl;
#endif
  int max_search= 350;

  for(int bi=startbin;bi<nbins;bi++){
    ANNpointArray       data_pts;        // data points to look for neighbors
    ANNpointArray       data_pts1;       // data points to measure dists for
    //int npts=ginbin[bi];
    int npts=pinbin[bi];
    if (npts == 0)
      continue;
    //cout<<"npts:"<<npts<<endl;
    int nsearch=ginbin[bi];
    //cout<<" bin has "<<npts<<" bright gals and "<<nsearch<<" total gals."<<endl;
    int k_search = 70;                   //initial search radius
    //if(npts<k_search) k_search=npts-1;
    if(nsearch<k_search) k_search=nsearch-1;
    nn_idx = new ANNidx[k_search];          // allocate near neigh indices
    dists = new ANNdist[k_search];          // allocate near neighbor dists

    data_pts1 = annAllocPts(npts, dim);         // allocate data points
    vector <int> ids(npts);
    vector <int> sids;
    sids.reserve(nsearch);
    if(bi>0) nsearch+=ginbin[bi-1];
    if(bi<ginbin.size()) nsearch+=ginbin[bi+1];
    data_pts = annAllocPts(nsearch, dim);         // allocate data points
    int nsi = 0;
    int npi = 0;

    for(int i=0;i<galaxies.size();i++){
      if((galaxies[i]->Zbin()==bi)||
	 (galaxies[i]->Zbin()==bi-1)||
	 (galaxies[i]->Zbin()==bi+1)){
	assert(nsi<nsearch);
	data_pts[nsi][0] = galaxies[i]->Ra();
	data_pts[nsi][1] = galaxies[i]->Dec();
	sids.push_back(i);
	nsi++;
      }
    }
    for(int i=0;i<points.size();i++){
      if(points[i]->Zbin()==bi){
	data_pts1[npi][0] = points[i]->Ra();
	data_pts1[npi][1] = points[i]->Dec();
	ids[npi]=i;
	npi++;
      } 
    }

    the_tree = new ANNkd_tree(
			      data_pts,  // the data points
			      nsearch,   // number of points
			      dim);		      // dimension of space

    for(int i=0;i<npts;i++){

      int gid = ids[i];
      #ifdef DEBUG
      if(pi%50000==0) {cout<<bi<<" "<<i<<" "<<pi<<endl; system("date");}
      #endif
      bool distfound = false;
      ANNpoint query_pt = data_pts1[i];
      float dist=0;
      int n_found = -1; 
      int t_found = -1; 
      while(!distfound){
	if(nsearch<k_search) k_search=nsearch-1;
	if(k_search>max_search) k_search=max_search;
	int nn = 0;

	the_tree->annkSearch(			// search
			     query_pt,		// query point
			     k_search,		// number of near neighbors
			     nn_idx,		// nearest neighbors (returned)
			     dists,		// distance (returned)
			     eps);		// error bound  
	for(int j=0;j<k_search;j++){
	  int chosen_id = sids[nn_idx[j]];

	  if((fabs(points[gid]->Z()-galaxies[chosen_id]->Z())<delta_z)){
	    ++nn;
	  }
	  if(nn==k_want){	

	    dist = points[gid]->ComovDist(sqrt(dists[j]));
	    n_found = j;
	    t_found = 0;
	    distfound = true;
	    break;
	  }
	}
	if(!distfound){

	  dist = points[pi]->ComovDist(sqrt(dists[k_search-1]));

	  if((k_search==nsearch-1)||(k_search>=max_search)){
	    n_found = k_search-1;
	    t_found = 1;
	    distfound = true;
	  }
	  if(dist>10){
	    n_found = -1;
	    t_found = 2;
	    distfound = true;
	  }
	  else{

	    nfout<<dist<<" "<<i<<" "<<points[i]->Z()<<" "<<k_search<<endl;

	    k_search = 2*k_search;
	    delete [] nn_idx;
	    delete [] dists;
	    nn_idx = new ANNidx[k_search];
	    dists = new ANNdist[k_search];
	  }
  	}//if dist is not found
      }//while dist is not found
      nndist[gid] = dist;
      ddout<<i<<" "<<pi<<" "<<gid<<" "<<dist<<" "<<t_found<<" "<<n_found<<" "<<k_search<<" "<<npts<<" "<<bi<<endl;
      pi++;
    }//loop over points in z bin
    //    ofstream ddout("dd2.out");
    //nnout<<nndist[i]<<endl;
    delete the_tree;
    delete [] nn_idx;
    delete [] dists;  
  }//loop over bins
  //MSG("Exiting neighbor");

  return nndist;
}



//vector <int> GetSEDs(vector <Galaxy *> &galaxies, vector <float> &nndist, vector <GalSED> & galseds){
vector <int> GetSEDs(vector <Galaxy *> &galaxies, vector <float> &nndist, vector <GalSED> & galseds, vector <Halo *> &halos){
  srand48(seed);
  vector <int> sed_ids(galaxies.size());
  float Percent = 0.0;
  double mr;
  cout<<galaxies[0]->Mr()<<endl;
  cout<<nndist[0]<<endl;

  //  omp_set_dynamic(0);
  //#pragma omp parallel for num_threads(4)
  for(int gi=0;gi<galaxies.size();gi++){
    //cout<<gi<<" "<<galaxies[gi]->Mr()<<endl;
    if(float(gi)/float(galaxies.size()) > Percent){
      //      cout<<"threads="<<omp_get_num_threads()<<endl;        
      cout<<"  "<<Percent*100.<<"% done"<<endl;
      Percent += 0.1;
    }

    Particle * p = galaxies[gi]->P();
    assert(p);

#ifdef DEBUG
    if(gi%20000==0) {cout<<gi<<endl; system("date");}
#endif
    mr = galaxies[gi]->Mr();
    sed_ids[gi] = findCloseGalaxies2(galseds, mr, nndist[gi], galaxies[gi]->Z(), galaxies[gi]->Central());
    if (gi==0)
      {
	cout << "Mr, dens found: " << galseds[sed_ids[0]].MR() << " " << galseds[sed_ids[0]].Dens() << endl;
      }
  }
  return sed_ids;
}

//This is for use without neighbor distances
//vector <int> GetSEDs(vector <Galaxy *> &galaxies, vector <GalSED> & galseds){
//  vector <int> sed_ids(galaxies.size());
//  for(int gi=0;gi<galaxies.size();gi++){
//    Particle * p = galaxies[gi]->P();
//    assert(p);  
//
//#ifdef DEBUG
//    if(gi%20000==0) {cout<<gi<<endl; system("date");}
//#endif
//    double mr = galaxies[gi]->Mr();
//    //if(mr>Magmin_col) mr = Magmin_col;
//
//    sed_ids[gi] = ChooseSED(galseds,mr,0, galaxies[gi]->Z(), galaxies[gi]->Central());
//    //}
//    //else sed_ids[gi]=-1;
//  }
//  return sed_ids;
//}

//This is for use without neighbor distances
vector <int> GetBCG_SEDs(vector <Galaxy *> &galaxies, vector <GalSED> & galseds){
  vector <int> sed_ids(galaxies.size());
  for(int gi=0;gi<galaxies.size();gi++){

      //}
    //else sed_ids[gi]=-1;
  }
  return sed_ids;
}

