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
  gadget_header head;
  long long npart;

  string infile = simulationfile;
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

  //binifstream pfile(infile.c_str());
  //  pfile>>buf;
  ifstream file (infile.c_str(), ios::in|ios::binary);
  if (file.is_open())
    {
      file.read((char *)(&buf), sizeof(buf));
      cout<<"buf: "<<buf<<endl;
      file.read((char *)(&head), sizeof(head));
      file.close();
      npart = head.npartTotal[1] + (((long long) head.npartTotalHighWord[0])<<32);
      np13 = int(pow((double) npart, ( double )1/3));

      Cosmology cosmo(float(head.Omega0), 0.0, float(head.HubbleParam));
      Simulation sim(simtype, cosmo);
      sim.Boxsize(head.BoxSize);
      sim.ParticleMass(head.mass[1]);
      sim.Np(np13);
    }
  if (buf != 256){
    cout<<"[DefineGadgetSimulation]  Incorrect header format for file "<<infile<<", buf = "<<buf<<endl;
    exit(1002);
  }
  //pfile.read((char *)(&head), sizeof(struct gadget_header));
  //  pfile.read((char *)(&head), sizeof(head));
  //  pfile.close();
  //  npart = head.npartTotal[1] + (((long long) head.npartTotalHighWord[0])<<32);
  //  np13 = int(pow((double) npart, ( double )1/3));

  //  Cosmology cosmo(float(head.Omega0), 0.0, float(head.HubbleParam));
  //  Simulation sim(simtype, cosmo);
  //  sim.Boxsize(head.BoxSize);
  //  sim.ParticleMass(head.mass[1]);
  //  sim.Np(np13);

  return sim;
}

Simulation DefineSimulation(void){
  Simulation simulation;
  if (simtype == "GADGET2") simulation=DefineGadgetSimulation();

  return simulation;
}


