#ifndef box_h
#define box_h
#include "simulation.h"
#include <string>
#include "ReadParameters.h"

//#define HVL
//#define LANL
#define GADGET
//#define MS_Gas
//#define CarmenLC

#ifdef JUST_COLORS
#ifndef COLORS
#define COLORS
#endif
//const std::string galaxy_file = "/nfs/slac/g/ki/ki01/mbusha/data/Alhambra/mock/code/junk_00L120_07.23.09.dat";
const std::string galaxy_file = "./galaxies.ascii";
#endif


//const shape_type shape = ZLIM;


static const int use_cells = 4; //number of square cells to select



static const float BOXFR = 1.0;

//in units of the boxx
//static const float XMIN = 0.0;
static const float XMIN = -1*BOXFR;
static const float XMAX = BOXFR;
//static const float YMIN = 0.0;
static const float YMIN = -1*BOXFR;
static const float YMAX = BOXFR;
//static const float ZMIN = 0.0;
static const float ZMIN = -1*BOXFR;
static const float ZMAX = BOXFR;
//note: for 5deg slice in dec centered on origin, 0.4927<z<0.5073.
//in Mpc/h
static const float RMIN = 0.0;
static const float RMAX = 3000;//1500.0;

//** careful.  you can't make this non-cubical without fixing number density stuff
//also if you want a sigma_8 shift much have xstart=ystart=zstart

#endif
