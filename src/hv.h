#ifndef hv_h
#define hv_h
#include "constants.h"
#include "particle.h"
#include "halo.h"
#include "stl_util.h" //for copy_if
#include "cluster.h"
#include "functions.h"

vector <Halo*> ReadHalos(void);
vector <Cluster*> ReadClusters(void);
//vector <Particle*> ReadParticles(int &nread);
vector <Particle*> ReadParticles();
vector <Particle*> ReadParticles(int &nread, vector <Halo *> halos);
//vector <Particle *> ReadWarren(void);

void HaloOcc(vector <Halo *> &halos);
#endif
