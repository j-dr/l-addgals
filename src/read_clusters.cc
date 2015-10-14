#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <list>
#include <strstream>

#include "biniostream.h"
#include "point.h"
#include "particle.h"
#include "cluster.h"
#include "singleton.h"
#include "functions.h"

vector <Cluster*> ReadClusters(){
  vector <Cluster*> clusters;

  string filename;
  filename = "/data1/risa/lens/"+simlabel+"_"+flabel+"_bcgs.dat";
  ifstream file(filename.c_str());
  if (file.fail()) {
    cerr<<"[read_centers error]: cannot open "<<filename<<endl;    exit(1);
  }
  else{
    cout<<"reading file "<<filename<<endl;
  }
  int hid = 0;
  while(file){
    float tt, r200, xx, yy, zz, vx, vy, vz, ra, dec, z;
    file>>tt>>xx>>yy>>zz>>vx>>vy>>vz>>ra>>dec>>z>>r200>>hid;
    int ngals = (int) (tt);
    //cout<<ngals<<" "<<xx<<" "<<ra<<" "<<dec<<" "<<hid<<"...";
    Point pos(xx, yy, zz);
    Point vel(vx, vy, vz);

    Cluster * cluster = new Cluster(pos, vel, ngals, ra, dec, z, r200, hid);
    //cout<<" "<<cluster->Ra()<<" "<<cluster->Dec()<<endl;
    clusters.push_back(cluster);
    hid++;
  }
  file.close();
  cout<<"[read_clusters] Read "<<clusters.size()<<" clusters."<<endl;
  return clusters;
}

