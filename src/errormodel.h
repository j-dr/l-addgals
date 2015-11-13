#ifndef ERRORMODEL_H
#define ERRORMODEL_H

#include <vector>

void apply_uniform_errormodel(float exptime[], float limmags[], float lnscat, int nband,
                              float zeropoint[], float nsigma, std::vector<float> &mag,
                              std::vector<float> &flux, std::vector<float> &fluxerr,
                              std::vector<float> &omag, std::vector<float> &omagerr);


void add_des_photometric_errors(float maglim[], int nband, std::vector<float> &mag,
				std::vector<float> &flux, std::vector<float> &fluxerr, 
				std::vector<float> &omag, std::vector<float> &omagerr,
				std::vector<bool> &good);

#endif
