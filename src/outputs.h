#ifndef OUTPUTS_H
#define OUTPUTS_H
#include <CCfits/CCfits>
#include <vector>
#include <galaxy.h>

void  write_bcc_catalogs(std::vector<Galaxy *> &galaxies, std::vector<Particle *> &particles,
			 std::vector<float> &amag, std::vector<float> &tmag, std::vector<float> &mr,
			 std::vector<float> &omag, std::vector<float> &omagerr, std::vector<float> &deltam,
			 std::vector<float> &flux, std::vector<float> &fluxerr, std::vector<double> &e,
			 std::vector<double> &s, std::vector<bool> &idx, std::vector <Halo *> &halos,
			 std::vector<int> &sed_ids, std::vector<float> &coeffs,
			 string outgfn, string outghfn);

void  write_bcc_catalogs_w_densities(std::vector<Galaxy *> &galaxies, std::vector<Particle *> &particles,
				     std::vector<float> &amag, std::vector<float> &tmag, std::vector<float> &mr,
				     std::vector<float> &omag, std::vector<float> &omagerr,
				     std::vector<float> &flux, std::vector<float> &deltam,
				     std::vector<float> &fluxerr, std::vector<double> &e,
				     std::vector<double> &s, std::vector<bool> &idx,std::vector <Halo *> &halos,
				     std::vector<int> &sed_ids, std::vector<float> &coeffs,
				     std::vector<float> &dist8,
				     std::vector<float> &nndist,
				     std::vector<float> &nndist_percent,
				     string outgfn, string outghfn);


#endif
