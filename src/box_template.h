#ifndef box_h
#define box_h
#include "simulation.h"
#include <string>

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

#ifdef HVL
const sim_type simulation = HV;
const hv_label box = HVLIGHT;
//const hv_label box = HVPOW;
const std::string simlabel = "hv";
const std::string simnum = "006";
#endif

#ifdef LANL
const sim_type simulation = WAR;
const war_label box = W384a;
const std::string simlabel = "wa";
const std::string simnum = "006";
#endif

#ifdef GADGET
const sim_type simulation = GADGET2;
//const gadget_label box = Millennium;
//const std::string simlabel = "downsample1000";
//const std::string simnum = "063";
const gadget_label box = Carmen;
//const std::string simlabel = "Carmen02_4004.le";
//const std::string simnum = "099";
const std::string simlabel = "PO";
//const std::string simlabel = "PO_downsample";
const std::string simnum = "000";
#endif

#ifdef CarmenLC
const sim_type simulation = GADGET2;
const gadget_label box = Carmen;
const std::string simlabel = "PO_wedge.2";
const std::string simnum = "000";
#endif

#ifdef MS_Gas
/*
const sim_type simulation = MGS;
const mgs_label box = LC;
//const gadget_label box = LC;
const std::string simlabel = "downsample100";
const std::string simnum = "063";
*/
const sim_type simulation = GADGET2;
const gadget_label box = Millennium;
const std::string simlabel = "MGS_lc";
const std::string simnum = "000";
#endif

enum shape_type{CUBE, ZLIM};
const shape_type shape = ZLIM;

#ifdef HVL
//for DES
//const std::string flabel = "POv1.08_test.01";//"v1.01_sm";//pow16na";//"pow16bbb";//v0.98n4_highz";//"v0.99c_2c_nntest1";//vhighz";//pdf";//_highz";
const std::string flabel = "MS_FullSky_v1.08_888.06";//MS_FullSky_788";
//const std::string flabel = "PO_dim_test";
const std::string datadir = "/nfs/slac/g/ki/ki01/mbusha/data/";
const std::string out_path = "/nfs/slac/g/ki/ki01/mbusha/data/";
const std::string path = out_path+"sdss/galcats/"+simlabel+"_"+flabel+"/";
#endif

#ifdef GADGET
//for ubercomp
//const std::string flabel = "v2.07";
//const std::string datadir = "/nfs/slac/g/ki/ki01/mbusha/projects/ubercomp/BaseData/Simulation/";
//const std::string out_path = "/nfs/slac/g/ki/ki01/mbusha/projects/ubercomp/BaseData/Simulation/mock/";
//const std::string path = out_path+simlabel+"_"+flabel+"/";

//for Millennium
//const std::string flabel = "sam_test_no_BCGs";
//const std::string datadir = "/nfs/slac/g/ki/ki02/mbusha/projects/addgals/Millennium/downsample1000/";
//const std::string out_path = "/nfs/slac/g/ki/ki01/mbusha/data/";
//const std::string path = out_path+"sdss/galcats/"+simlabel+"_"+flabel+"/";

//for LasDamas
//const std::string flabel = "Box120";
///const std::string flabel = "Carmen02_box_999";
//const std::string datadir = "/nfs/slac/g/ki/ki04/mbusha/projects/LasDamas/Simulations/Carmen/02/";
const std::string flabel = "Carmen02_999";
const std::string datadir = "/nfs/slac/g/ki/ki04/mbusha/projects/LasDamas/Simulations/Carmen/02/analysis/LC/PO_no_tile/";
//const std::string datadir = "/nfs/slac/g/ki/ki04/mbusha/projects/LasDamas/Simulations/Carmen/02/analysis/LC/PO_no_tile/downsample/";
const std::string out_path = "./hv_output/";
const std::string path = out_path;
#endif

#ifdef CarmenLC
const std::string flabel = "Carmen_dc5_100";
const std::string datadir = "/nfs/slac/g/ki/ki04/mbusha/projects/LasDamas/Simulations/Carmen/02/analysis/LC/dc5/";
//const std::string out_path = "/nfs/slac/g/ki/ki01/mbusha/data/";
//const std::string path = out_path+"sdss/galcats/"+simlabel+"_"+flabel+"/";
const std::string out_path = "./hv_output/";
const std::string path = out_path;
#endif

#ifdef MS_Gas
const std::string flabel = "ADDGALS_MGS";
const std::string datadir = "/nfs/slac/g/ki/ki02/mbusha/projects/addgals/Millennium/mgs/data/mgs_lc_gadget/";
const std::string out_path = "./hv_output/";
const std::string path = out_path;
#endif

////path for output files
//const std::string datadir = "/data/risa/";
//const std::string datadir = "/usr/work/risa/";
//const std::string datadir = "/home/risa/data/";

static const int use_cells = 4; //number of square cells to select

//angle limits
//note: dec ranges from -90 to 90.  ra ranges from -180 to 180.//0 to 360.
//positive x corresponds to -90<ra<90.
//positive y corresponds to 0<ra<180.
//positive z corresponds to postive dec9

//static const float DECMIN = 35.;
//static const float DECMAX = 50.;//16.0;
//static const float RAMIN = 10.;
//static const float RAMAX = 30.0;//32.0;
static const float DECMIN = 0.0;
static const float DECMAX = 90.0;
static const float RAMIN = 0.0;
static const float RAMAX = 90.0;

//redshift limit
#ifdef PARALLEL
extern float ZREDMIN;
extern float ZREDMAX;
static const float SIM_ZREDMIN = 0.;
static const float SIM_ZREDMAX = 0.3;
extern int NTasks;
extern int ThisTask;
#else
static const float ZREDMIN
//static const float ZREDMAX = 0.10; //small cube
//static const float ZREDMAX = 0.064; //one cube
//static const float ZREDMAX = 0.129; //two cubes
//static const float ZREDMAX = 0.196; //three cubes
//static const float ZREDMAX = 0.266; //four cubes
//static const float ZREDMAX = 0.34;  //five cubes
//static const float ZREDMAX = 0.416; //six cubes
//static const float ZREDMAX = 0.496; //seven cubes
//static const float ZREDMAX = 0.57; //eight cubes
//static const float ZREDMAX = 1.4;
//static const float ZREDMAX = 3.0;
static const float ZREDMAX
//static const float ZREDMAX = 0.128; //ubercomp
//static const float ZREDMAX = 0.31;
#endif

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
