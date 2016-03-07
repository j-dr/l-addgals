#ifndef OUTPUTS_H
#define OUTPUTS_H
#include <CCfits/CCfits>
#include <vector>
#include <galaxy.h>

void write_bcc_catalogs(std::vector<float> &amag, std::vector<float> &tmag, 
			std::vector<float> &sdssr, std::vector<float> &omag, 
			std::vector<float> &omagerr, std::vector<float> &flux, 
			std::vector<float> &ivar, std::vector<double> &e, std::vector<double> &s,
			std::vector<float> &ra, std::vector<float> &dec, std::vector<float> &px,
			std::vector<float> &py, std::vector<float> &pz, std::vector<float> &vx,
			std::vector<float> &vy, std::vector<float> &vz, std::vector<float> &z,
			std::vector<int> &gid, std::vector<int> &central, 
			std::vector<int> &haloid, std::vector<int> &ecatid, 
			std::vector<float> &coeffs, std::string outgfn, std::string outghfn);

int extract_galaxy_information(std::vector<Galaxy *> &galaxies, std::vector<Particle *> &particles, 
			       std::vector<float> &amag, std::vector<float> &mr, 
			       std::vector<bool> &idx, std::vector <Halo *> &halos, 
			       std::vector<float> &coeffs, std::vector<float> &ra, 
			       std::vector<float> &dec, std::vector<float> &px, 
			       std::vector<float> &py, std::vector<float> &pz, 
			       std::vector<float> &vx, std::vector<float> &vy, 
			       std::vector<float> &vz, std::vector<float> &sdssr, 
			       std::vector<float> &z, std::vector<int> &gid, 
			       std::vector<int> &central, std::vector<int> &haloid, 
			       std::vector<float> &cf, std::vector<float> &am, 
			       std::vector<int> &ecatid, int nbands, int ntemp);
#endif









