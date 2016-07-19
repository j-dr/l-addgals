#ifndef GALAXY_GLOBAL_H
#define GALAXY_GLOBAL_H
#include <vector>

#define NBIN 8500
struct den_ent{
  float prob[NBIN];
  float r[NBIN];
};

void read_out_galaxy_info(vector<Galaxy *> &gal, vector<float> &mr,
			  vector<float> &z, vector<GalSED> &seds,
			  vector<int> &sed_id, vector<int> &sed_cat_id);
void read_out_galaxy_info_w_densities(vector<Galaxy *> &gal,
              vector<float> &mr,
              vector<float> &z, vector<GalSED> &seds,
              vector<int> &sed_id, vector<int> &sed_cat_id,
              vector<float> &dist8);
float LocalDens(den_ent pdf);
float SelectGalaxyZ();
std::vector <Galaxy *> GetGalaxies(double vol, float phi_rescale);

#endif
