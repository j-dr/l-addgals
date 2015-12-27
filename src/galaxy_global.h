#ifndef GALAXY_GLOBAL_H
#define GALAXY_GLOBAL_H
#include <vector>

#define NBIN 8500
struct den_ent{
  float prob[NBIN];
  float r[NBIN];
};

void read_out_galaxy_info(vector<Galaxy *> &gal, vector<GalSED> &sed,
                          vector<int> &sed_id, vector<float> &mr,
                          vector<float> &z, vector<int> &id);
float LocalDens(den_ent pdf);
float SelectGalaxyZ();
std::vector <Galaxy *> GetGalaxies(double vol, float phi_rescale);

#endif
