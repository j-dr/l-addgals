#include "color.h"

vector <GalSED> ReadSED(){
  vector <GalSED> galseds;
  galseds.reserve(38358);
  //   string filename = "../sdssinp/dr3sedv3.2.dat";
#ifndef COLORS_FROM_RELATIVE_DENSITY
   string filename = "../sdssinp/dr3_vagc_id.dat";
#else
   /*
#ifdef RED_FRACTION
   string filename = "/nfs/slac/g/ki/ki01/mbusha/data/sdss/dr6/cooper/dr6_cooper_id_with_red.dat";
#else
   string filename = "/nfs/slac/g/ki/ki01/mbusha/data/sdss/dr6/cooper/dr6_cooper_id.dat";
#endif
   */
   string filename = colortrdir + "/dr6_cooper_id_with_red.dat";
#endif

   ifstream file(filename.c_str());
   if (file.fail()) {
     cerr<<"error: cannot open "<<filename<<endl;
     exit(1);
   }
   int badgals=0;
   int goodgals=0;

   int catid = 0;
   while(file){
     //long int id;
     double Mr, d10, d5, dd;
     int red;
     file>>Mr>>d10>>d5>>red;
#ifndef RED_FRACTION
     red = 0;
#endif
      catid++;
    if((Mr<=Magmin)&&(d10>0)){
      if(dens_measure == TENTH)
	dd = d10;
      else
	dd = d5;
      GalSED galsed(Mr,dd,goodgals,catid,red);
      /*
#ifdef RED_FRACTION
      if (catid < 20){
	cout<<galsed.MR()<<" "<<galsed.Dens()<<" "<<galsed.Red()<<" "<<galsed.Id()<<endl;
      }
#endif
      */
      //galsed.Id(goodgals); //just the number that have been used so far
      galseds.push_back(galsed);
      //galseds.CatId(catid);
      //      SEDTuple SED(coef[0],coef[1],coef[2],Mr, sdssid);
      //      SEDs.push_back(SED);
      goodgals++;
      d10 = 0;
    }
    else {
      //      cout<<Magmin<<" "<<Mr<<" "<<" "<<d10<<" "<<d5<<endl;
      badgals++;
    }
  }
  file.close();
  //  cout<<SEDs.size()<<" ";
  //SEDs[SEDs.size()-1].Print();
  return galseds;
}


vector <GalSED> ReadDimSED(){
  vector <GalSED> galseds;
  galseds.reserve(150000);
  //   string filename = "../sdssinp/dr3sedv3.2.dat";
  string filename = "../sdssinp/dr3_dim_vagc_id.dat";

  ifstream file(filename.c_str());
  if (file.fail()) {
    cerr<<"error: cannot open "<<filename<<endl;
    exit(1);
  }
  int badgals=0;
  int goodgals=0;
  
  long int catid = 0;
  while(file){
    double Mr;
    file>>Mr;
    catid++;
    GalSED galsed(Mr,0,goodgals,catid);
    galseds.push_back(galsed);
    goodgals++;
  }
  file.close();
  return galseds;
}

#ifdef DR3_SED
vector <GalSED> ReadSED(vector <SEDTuple> &SEDs){
  vector <GalSED> galseds;
  galseds.reserve(28881);
  //   string filename = "../sdssinp/dr3sedv3.2.dat";
   string filename = "../sdssinp/dr3_vagc.dat";

   ifstream file(filename.c_str());
   if (file.fail()) {
     cerr<<"error: cannot open "<<filename<<endl;
     exit(1);
   }
   int badgals=0;
   int goodgals=0;

   while(file){
     vector <double> coef(3);
     vector <int> sdssid(3);
     double Mr, d10, d5, dd;// plate, mjd, fiber;
     //file>>coef[0]>>coef[1]>>coef[2]>>coef[3]>>Mr>>d10>>d5;
    //sdss id is currently plate mjd, fiber
     //    file>>Mr>>coef[0]>>coef[1]>>coef[2]>>d10>>d5>>sdssid[0]>>sdssid[1]>>sdssid[2];
     file>>id>>Mr>>d10>>d5;
    if((Mr<=Magmin)&&(d10>0)){
      if(dens_measure == TENTH)
	dd = d10;
      else
	dd = d5;
      GalSED galsed(Mr,dd);
      galsed.Id(goodgals); //just the number that have been used so far
      galseds.push_back(galsed);
      SEDTuple SED(coef[0],coef[1],coef[2],Mr, sdssid);
      SEDs.push_back(SED);
      goodgals++;
      d10 = 0;
    }
    else {
      //      cout<<Magmin<<" "<<Mr<<" "<<" "<<d10<<" "<<d5<<endl;
      badgals++;
    }
  }
  file.close();
  cout<<SEDs.size()<<" ";
  SEDs[SEDs.size()-1].Print();
  return galseds;
}
#endif

#ifdef DR2
vector <GalSED> ReadSED(vector <SEDTuple> &SEDs){
  vector <GalSED> galseds;
  galseds.reserve(19666);
  string filename = "../sdssinp/dr2bsedv3.2.dat";
  //string filename = "../sdssinp/dr1sedv3.2.dat";
  //string filename = "../sdssinp/dr1sed.dat";
  ifstream file(filename.c_str());
  if (file.fail()) {
    cerr<<"error: cannot open "<<filename<<endl;
    exit(1);
  }
  int badgals=0;
  int goodgals=0;

  while(file){
    vector <double> coef(3);
    double Mr, d10, d5, dd;
    //file>>coef[0]>>coef[1]>>coef[2]>>coef[3]>>Mr>>d10>>d5;
    file>>Mr>>coef[0]>>coef[1]>>coef[2]>>d10>>d5;
    if((Mr<=Magmin)){
      if(dens_measure == TENTH)
	dd = d10;
      else
	dd = d5;
      GalSED galsed(Mr,dd);
      galsed.Id(goodgals); //just the number that have been used so far
      galseds.push_back(galsed);
      SEDTuple SED(coef[0],coef[1],coef[2],Mr);
      SEDs.push_back(SED);
      goodgals++;
    }
    else {
      //      cout<<Magmin<<" "<<Mr<<" "<<" "<<d10<<" "<<d5<<endl;
      badgals++;
    }
  }
  file.close();
  cout<<"rejected "<<badgals<<" dim galaxies, using "<<goodgals<<" "<<galseds.size()<<endl;
  return galseds;
}

#endif
