#include <iostream>
#include <fstream>
#include <vector>
#include <ANN/ANN.h>			// ANN declarations
#include "nr.h"
#include "galaxy.h"
#include "singleton.h"
#include "choose.h"
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
  cout<<"spherical gnp"<<endl;
  vector <float> color_percent(galaxies.size());
  int nColorBins = floor(zmax->GetVal()/ColorBinSize)+1;
  struct GalaxyPercent tGalaxyId(0.,0);
  vector <GalaxyPercent> GalaxyId;
  cout<<"setup donem now looping through galaxies"<<endl;
  for (int i=0;i<galaxies.size();i++)
    {
      float tnndist = nndist[i];
      tGalaxyId.dens(tnndist);
      tGalaxyId.gid(i);
      GalaxyId.push_back(tGalaxyId);
    }
  cout<<"now sorting"<<endl;
  sort(GalaxyId.begin(),GalaxyId.end(),CompareByDens);
  cout<<"calculating percents"<<endl;
  for (int i=0;i<galaxies.size();i++) 
    color_percent[GalaxyId[i].gid()] = float(i)/float(galaxies.size());
  cout<<"clearing galaxyid"<<endl;
  GalaxyId.clear();

  float min_percent = 100.0;
  float max_percent = 0.0;
  for (int i=0;i<galaxies.size();i++){
    if (color_percent[i] < min_percent) min_percent = color_percent[i];
    if (color_percent[i] > max_percent) max_percent = color_percent[i];
  }

  cout<<"First 3 color percents: "<<color_percent[0]<<" "<<color_percent[1]<<" "<<color_percent[2]<<endl;
  cout<<"Last 3 color percents: "<<color_percent[galaxies.size()-1]<<" "<<color_percent[galaxies.size()-2]<<" "<<color_percent[galaxies.size()-1]<<endl;
  cout<<"min/max percent: "<<min_percent<<" "<<max_percent<<endl;

  return color_percent;
}
#endif

void printPt(ANNpoint p)
{
  cout<<"("<<p[0];
  for (int i=1;i<3;i++) cout<<", "<<p[i];
  cout<<")\n";
}

//calculate the number of galaxies within 10 Mpc/h
vector <float> GetNeighborDist(vector <Galaxy *> Junk, vector <Galaxy *> galaxies){
  cout<<"In GetNeighborNum"<<endl;

  int maxPts = galaxies.size(); //actual number of data points
  int nPts;
  ANNpointArray dataPts; //array of data points
  ANNpoint queryPt; //query point
  ANNkd_tree * kdTree;
  ANNdist r_search = 10.0; //our search radius
  ANNdist r2_search = r_search*r_search;
  ANNdist r3_search = r2_search*r_search;
  int knn = 0; //max number of neighbors in search radius
  double eps = 0.0;
  int dim = 3;
  vector <float> nndist(galaxies.size());


  //allocate our arrays for the galaxy positions
  queryPt = annAllocPt(dim);
  dataPts = annAllocPts(maxPts,dim);
  nPts = 0;
  for(int ig=0;ig<maxPts;ig++)
    {
      if (galaxies[ig]->Mr() > Magmin_dens) continue;
      dataPts[nPts][0] = galaxies[ig]->P()->X();
      dataPts[nPts][1] = galaxies[ig]->P()->Y();
      dataPts[nPts][2] = galaxies[ig]->P()->Z();
      if (nPts < 10) {printPt(dataPts[nPts]);}
      nPts++;
    }
  cout<<"Number of galaxies brighter than "<<Magmin_dens<<": "<<nPts<<endl;

  //build our tree
  kdTree = new ANNkd_tree(dataPts, nPts, dim);

  //find densities of each galaxy
  cout<<"Finding the densitities..."<<endl;
  for(int ig=0;ig<galaxies.size();ig++)
    {
      //copy position to query point
      //cout<<"Filling the query pt..."<<endl;
      queryPt[0] = galaxies[ig]->P()->X();
      queryPt[1] = galaxies[ig]->P()->Y();
      queryPt[2] = galaxies[ig]->P()->Z();

      //do the query
      nndist[ig] = kdTree->annkFRSearch(queryPt, r2_search, knn);

      if (ig < 10) cout<<r2_search<<" "<<nndist[ig]<<endl;
      
    }
  int min_nndist = 1000;
  int max_nndist = 0;
  int nzero = 0;
  for (int ig=0;ig<galaxies.size();ig++){
    if (nndist[ig] < min_nndist) min_nndist = nndist[ig];
    if (nndist[ig] > max_nndist) max_nndist = nndist[ig];
    if (nndist[ig] == 0) {nzero++; galaxies[ig]->P()->PosPrint();}
  }
  cout<<"min/max measured number of objects within 10 Mpc: "<<min_nndist<<" "<<max_nndist<<endl;

  cout<<"Number of objects with 0 neighbors: "<<nzero<<endl;

  delete kdTree;
  annClose();

  return nndist;
}



//vector <int> GetSEDs(vector <Galaxy *> &galaxies, vector <float> &nndist, vector <GalSED> & galseds){
vector <int> GetSEDs(vector <Galaxy *> &galaxies, vector <float> &nndist, vector <GalSED> & galseds, vector <Halo *> &halos){
  srand48(seed);
  vector <int> sed_ids(galaxies.size());
  float Percent = 0.0;
  cout<<galaxies[0]->Mr()<<endl;

  for(int gi=0;gi<galaxies.size();gi++){
    //cout<<gi<<" "<<galaxies[gi]->Mr()<<endl;
    if(float(gi)/float(galaxies.size()) > Percent){
      cout<<"  "<<Percent*100.<<"% done"<<endl;
      Percent += 0.1;
    }

    Particle * p = galaxies[gi]->P();
    assert(p);  //this better be true since you removed the other ones.
#ifdef DEBUG
    if(gi%20000==0) {cout<<gi<<endl; system("date");}
#endif
    double mr = galaxies[gi]->Mr();

    sed_ids[gi] = findCloseGalaxies2(galseds, mr, nndist[gi], galaxies[gi]->Z(), galaxies[gi]->Central());

  }
  return sed_ids;
}

//This is for use without neighbor distances
vector <int> GetSEDs(vector <Galaxy *> &galaxies, vector <GalSED> & galseds){
  vector <int> sed_ids(galaxies.size());
  for(int gi=0;gi<galaxies.size();gi++){
    Particle * p = galaxies[gi]->P();
    assert(p);  //this better be true since you removed the other ones.
#ifdef DEBUG
    if(gi%20000==0) {cout<<gi<<endl; system("date");}
#endif
    double mr = galaxies[gi]->Mr();

    sed_ids[gi] = ChooseSED(galseds,mr,0, galaxies[gi]->Z(), galaxies[gi]->Central());
  }
  return sed_ids;
}

//This is for use without neighbor distances
vector <int> GetBCG_SEDs(vector <Galaxy *> &galaxies, vector <GalSED> & galseds){
  vector <int> sed_ids(galaxies.size());
  for(int gi=0;gi<galaxies.size();gi++){

      //}
    //else sed_ids[gi]=-1;
  }
  return sed_ids;
}

