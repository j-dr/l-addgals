#ifndef ERRORMODEL_H
#define ERRORMODEL_H

#include <vector>

void apply_uniform_errormodel(float exptime[], float limmags[], float lnscat[], 
			      int nband, float zeropoint[], float nsigma, 
			      std::vector<float> &mag, std::vector<float> &flux, 
			      std::vector<float> &fluxerr, std::vector<float> &omag, 
			      std::vector<float> &omagerr);

void observe_des_y5(std::vector<float> &mag, std::vector<float> &flux, 
		    std::vector<float> &fluxerr, std::vector<float> &omag,
		    std::vector<float> &omagerr);

#endif
