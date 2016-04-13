#include <string>
#include <iostream>
#include <cassert>
#include <fstream>
#include "biniostream.h"
#include "global_vars.h"
#include "ReadParameters.h"
#include "gadget.h"

Simulation DefineGadgetSimulation(void){

  int buf, np13;
  long long npart;
  string infile = simulationfile;
  gadget_header head;
  const char *infilechar = infile.c_str();
  ifstream inFile(infilechar, ios::in);
  if(!inFile){
    inFile.close();
    infile += ".0";
    const char *infilechar2 = infile.c_str();
    ifstream inFile(infilechar2, ios::in);
    if(!inFile){
      cout<<"Error!  The input file does not exist: "<<infile<<endl;
      exit(1001);
    }
  }
  inFile.close();

  binifstream pfile(infile.c_str());
  pfile>>buf;
  if (buf != 256){
    cout<<"Incorrect header format for file "<<infile<<", buf = "<<buf<<endl;
    exit(1002);
  }
  //pfile.read((char *)(&head), sizeof(struct gadget_header));
  pfile.read((char *)(&head), sizeof(head));
  pfile.close();
  npart = head.npartTotal[1] + (((long long) head.npartTotalHighWord[0])<<32);
  np13 = int(pow((double) npart, ( double )1/3));

  Cosmology cosmo(float(head.Omega0), 0.0, float(head.HubbleParam));
  Simulation sim(simtype, cosmo);
  sim.Boxsize(head.BoxSize);
  sim.ParticleMass(head.mass[1]);
  cout << "particle mass: " << sim.ParticleMass() << endl;
  sim.Np(np13);

  return sim;
}


Simulation DefineBCCGadgetSimulation(void){

  int buf, np13;
  long long npart;
  io_header head;
  int r = 0;

  while (true) {
    std::ostringstream convert;
    convert << datadir << simlabel << "_000_" << r << "_0";
    std::string fname = convert.str();
    std::ifstream pfile(fname.c_str());
    
    if (pfile.fail()) {
      r++;
      continue;
    }
    
    pfile.read((char *)(&head), sizeof(struct io_header));
    break;
  }

  Cosmology cosmo(head.Omega0, 0.0, head.HubbleParam);
  Simulation sim(simtype, cosmo);
  sim.Boxsize(head.BoxSize);
  sim.ParticleMass(head.mass*pow(10.0,10.0));
  sim.Np(0);

  return sim;
}



Simulation DefineSimulation(void){
  Simulation simulation;
  if (simtype == "GADGET2") simulation=DefineGadgetSimulation();
  if (simtype == "BCCGADGET2") simulation=DefineBCCGadgetSimulation();

  return simulation;
}


