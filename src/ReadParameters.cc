#include "ParameterDatabase.h"
#include "StringDatabase.h"

 string simtype; 	//Gadget or other data format
 string simlabel;	//Not currently being used
// string simnum ;	//Not currently bing used
// string flabel ;	//Not currently being used
 string datadir ;	//Only used in the (incomplete) MGS readin function
 string halofile ;	//File where halodata is sotred
 string rnn_halofile ;	//File where halodata is sotred
 string simulationfile; //file with particle data
 string rnnfile;
 string out_path ;
 string denspdffile ;
 string lbcgfile ;
 float DECMIN ;
 float DECMAX;
 float RAMIN ;
 float RAMAX;
 float ZREDMIN;
 float ZREDMAX ;
 float sinDECMAX;
 float sinDECMIN;
 float RMIN_REAL;
 float RMAX_REAL;
 float Z_REDFRACTION1;
 float Z_REDFRACTION2;
 float REDFRACTION1;
 float REDFRACTION2;
 float SCATTER;
 int REDSHIFT_FIT;
 int GLOBAL_FIT;
 float cm0, cm1, cm2, cm3, cm4, cmz1, cmz2, cmz3;
 float cs0, cs1, cs2, cs3, csz1, csz2;
 float fm0, fm1, fm2, fm3, fmz1, fmz2;
 float fs0, fs1, fs2, fs3, fs4, fsz1, fsz2, fsz3;
 float p0, p1, p2, p3, pz1, pz2, pz3;
#ifdef HEALPIX
 long nSide;
 long PixelNum;
#endif
 float minrnn, maxrnn;
#ifdef SHAM_TEST
 string sham_file;
#endif
#ifdef BCC
 int LCNUM;
#endif
 double sim_redshift;
  
  //hod file parameters
  string hodfile;
  float mhost_cut;
  float rnn_cut;
  int read_hod = 0;

 float ngals ; // Was int actually
 float PARTICLE_FRACTION;
 //enum ev_type{NOEV, BLAN, FABER, TIME};
 //ev_type evolution = BLAN; // Was enum actually
 ev_type evolution = TIME; // Was enum actually
 float Mstar;
 float phistar;
 float alpha;
 float Q;
 float Magmin;
 float Magmin_dens;
 float oMagmin ;
// float Magmin_pdf;
 float BCG_Mass_lim;
 float zTol; 
 float ColorBinSize;


void fillParameters(ParameterDatabase* pd)
{
	double angle_const = 45.0/atan(1.0);

	DECMIN = pd->findParameterValue("DECMIN");
	sinDECMIN = sin(DECMIN/angle_const);
	DECMAX = pd->findParameterValue("DECMAX");
	sinDECMAX = sin(DECMAX/angle_const);
	RAMIN = pd->findParameterValue("RAMIN");
	RAMAX = pd->findParameterValue("RAMAX");
	ZREDMIN = pd->findParameterValue("ZREDMIN");
	ZREDMAX = pd->findParameterValue("ZREDMAX");
        REDFRACTION1 = pd->findParameterValue("REDFRACTION1");
        REDFRACTION2 = pd->findParameterValue("REDFRACTION2");
        Z_REDFRACTION1 = pd->findParameterValue("Z_REDFRACTION1");
        Z_REDFRACTION2 = pd->findParameterValue("Z_REDFRACTION2");
        SCATTER = pd->findParameterValue("SCATTER");
	REDSHIFT_FIT = 0;
        //REDSHIFT_FIT = pd->findParameterValue("GlobalFit");
#ifdef BCC
	LCNUM = pd->findParameterValue("LCNUM");
#endif
	GLOBAL_FIT=1;
	if(GLOBAL_FIT){
		cm0 = pd->findParameterValue("cm0");
		cm1 = pd->findParameterValue("cm1");
		cm2 = pd->findParameterValue("cm2");
		cm3 = pd->findParameterValue("cm3");
		cm4 = pd->findParameterValue("cm4");
		cmz1 = pd->findParameterValue("cmz1");
		cmz2 = pd->findParameterValue("cmz2");
		cmz3 = pd->findParameterValue("cmz3");
		cs0 = pd->findParameterValue("cs0");
		cs1 = pd->findParameterValue("cs1");
		cs2 = pd->findParameterValue("cs2");
		cs3 = pd->findParameterValue("cs3");
		csz1 = pd->findParameterValue("csz1");
		csz2 = pd->findParameterValue("csz2");
		fm0 = pd->findParameterValue("fm0");
		fm1 = pd->findParameterValue("fm1");
		fm2 = pd->findParameterValue("fm2");
		fm3 = pd->findParameterValue("fm3");
		fmz1 = pd->findParameterValue("fmz1");
		fmz2 = pd->findParameterValue("fmz2");
		fs0 = pd->findParameterValue("fs0");
		fs1 = pd->findParameterValue("fs1");
		fs2 = pd->findParameterValue("fs2");
		fs3 = pd->findParameterValue("fs3");
		fs4 = pd->findParameterValue("fs4");
		fsz1 = pd->findParameterValue("fsz1");
		fsz2 = pd->findParameterValue("fsz2");
		fsz3 = pd->findParameterValue("fsz3");
		p0 = pd->findParameterValue("p0");
		p1 = pd->findParameterValue("p1");
		p2 = pd->findParameterValue("p2");
		p3 = pd->findParameterValue("p3");
		pz1 = pd->findParameterValue("pz1");
		pz2 = pd->findParameterValue("pz2");
		pz3 = pd->findParameterValue("pz3");
	}
#ifdef HEALPIX
	nSide = pd->findParameterValue("nSide");
	PixelNum = pd->findParameterValue("PixelNum");
#endif

	ngals = pd->findParameterValue("ngals");
	PARTICLE_FRACTION = pd->findParameterValue("PARTICLE_FRACTION");
	evolution = (ev_type) ((int)pd->findParameterValue("evolution"));
	Mstar = pd->findParameterValue("Mstar");
	phistar = pd->findParameterValue("phistar");
	alpha = pd->findParameterValue("alpha");
	Q = pd->findParameterValue("Q");
	Magmin = pd->findParameterValue("Magmin");
	Magmin_dens = pd->findParameterValue("Magmin_dens");
	oMagmin = pd->findParameterValue("oMagmin");
//	Magmin_pdf = pd->findParameterValue("Magmin_pdf");
	BCG_Mass_lim = pd->findParameterValue("BCG_Mass_lim");
	zTol = pd->findParameterValue("zTol"); 
	ColorBinSize = pd->findParameterValue("ColorBinSize");
	mhost_cut = pd->findParameterValue("mhost_cut");
	if (mhost_cut == -1) mhost_cut = 1e20;
	rnn_cut = pd->findParameterValue("rnn_cut");
}

void fillParameters(StringDatabase* sd)
{
	simtype = sd->findParameterValue("simtype");
	simlabel = sd->findParameterValue("simlabel");
	//simnum = sd->findParameterValue("simnum");
	//flabel = sd->findParameterValue("flabel");
	datadir = sd->findParameterValue("datadir");
	halofile = sd->findParameterValue("halofile");
	rnn_halofile = sd->findParameterValue("rnn_halofile");
	simulationfile = sd->findParameterValue("simulationfile");
	rnnfile = sd->findParameterValue("rnnfile");
	out_path = sd->findParameterValue("out_path");
	denspdffile = sd->findParameterValue("denspdffile");
	lbcgfile = sd->findParameterValue("lbcgfile");
#ifdef SHAM_TEST
	sham_file = sd->findParameterValue("SHAM_file");
#endif
	hodfile = sd->findParameterValue("hodfile");
	if (hodfile != "not_found") read_hod = 1;
}


void readParameters()
{
	//Give the name of file which contains numerical parameters only
	ParameterDatabase* pd = new ParameterDatabase("NumericalParameters");
	//Give the name of file which contains string parameters only
	StringDatabase* sd = new StringDatabase("StringParameters");
	fillParameters(pd);
	fillParameters(sd);
	//To see the list of parameters, uncomment the two statements below
	//pd->print();
	//sd->print();
}
