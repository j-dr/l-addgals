#ifndef read_h
#define read_h
#include <string>

vector <Halo*> ReadWHalos(void);
vector <Halo*> ReadHVHalos(void);
vector <Particle *> ReadWarren(void);
vector <Particle *> ReadHVParticles(void);

int readFitsHeader(std::string filename);

#endif
