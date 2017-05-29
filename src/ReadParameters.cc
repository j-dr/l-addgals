#include "ParameterDatabase.h"
#include "StringDatabase.h"

 string simtype; 	//Gadget or other data format
 string simlabel;	//Not currently being used
// string simnum ;	//Not currently bing used
 string flabel ;	
 string datadir ;	
 string halofile ;	//File where halodata is sotred
 string rnn_halofile ;	//File where halodata is sotred
 string simulationfile; //file with particle data
 string rnnfile;
 string out_path ;
 string denspdffile ;
 string lbcgfile ;
 string colortrdir;
 int EVOLVECEN;
 float QBASELINE;
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
 int LCNUM;
 float cm0, cm1, cm2, cm3, cm4, cmz1, cmz2, cmz3, cmz4, cm1z1, cm1z2, cm2z1, cm2z2, cm3z1, cm1z3;
 float cs0, cs1, cs2, cs3, cs4, csz1, csz2, csz3, csz4, cs1z1, cs1z2, cs2z1, cs2z2, cs3z1, cs1z3;
 float fm0, fm1, fm2, fm3, fm4, fmz1, fmz2, fmz3, fmz4, fm1z1, fm1z2, fm2z1, fm2z2, fm3z1, fm1z3;
 float fs0, fs1, fs2, fs3, fs4, fsz1, fsz2, fsz3, fsz4, fs1z1, fs1z2, fs2z1, fs2z2, fs3z1, fs1z3;
 float p0, p1, p2, p3, p4, pz1, pz2, pz3, pz4, p1z1, p1z2, p2z1, p2z2, p3z1, p1z3;
#ifdef HEALPIX
 long nSide;
 long PixelNum;
#endif
 float minrnn, maxrnn;
#ifdef SHAM_TEST
 string sham_file;
#endif
#ifdef BCC
 string PSTR;
 string ZSTR;
#endif
#ifdef FITS_GALAXIES
 string ffile;
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
	EVOLVECEN = pd->findParameterValue("EVOLVECEN");
	QBASELINE = pd->findParameterValue("QBASELINE");
	REDSHIFT_FIT = 0;
	if (QBASELINE==-1) QBASELINE = 0;
	if (EVOLVECEN==-1) EVOLVECEN = 1;
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
		cmz4 = pd->findParameterValue("cmz4");
		cm1z1 = pd->findParameterValue("cm1z1");
		cm1z2 = pd->findParameterValue("cm1z2");
		cm2z1 = pd->findParameterValue("cm2z1");
		cm2z2 = pd->findParameterValue("cm2z2");
		cm3z1 = pd->findParameterValue("cm3z1");
		cm1z3 = pd->findParameterValue("cm1z3");
		cs0 = pd->findParameterValue("cs0");
		cs1 = pd->findParameterValue("cs1");
		cs2 = pd->findParameterValue("cs2");
		cs3 = pd->findParameterValue("cs3");
		cs4 = pd->findParameterValue("cs4");
		csz1 = pd->findParameterValue("csz1");
		csz2 = pd->findParameterValue("csz2");
		csz3 = pd->findParameterValue("csz3");
		csz4 = pd->findParameterValue("csz4");
		cs1z1 = pd->findParameterValue("cs1z1");
		cs1z2 = pd->findParameterValue("cs1z2");
		cs2z1 = pd->findParameterValue("cs2z1");
		cs2z2 = pd->findParameterValue("cs2z2");
		cs3z1 = pd->findParameterValue("cs3z1");
		cs1z3 = pd->findParameterValue("cs1z3");
		fm0 = pd->findParameterValue("fm0");
		fm1 = pd->findParameterValue("fm1");
		fm2 = pd->findParameterValue("fm2");
		fm3 = pd->findParameterValue("fm3");
		fm4 = pd->findParameterValue("fm4");
		fmz1 = pd->findParameterValue("fmz1");
		fmz2 = pd->findParameterValue("fmz2");
		fmz3 = pd->findParameterValue("fmz3");
		fmz4 = pd->findParameterValue("fmz4");
		fm1z1 = pd->findParameterValue("fm1z1");
		fm1z2 = pd->findParameterValue("fm1z2");
		fm2z1 = pd->findParameterValue("fm2z1");
		fm2z2 = pd->findParameterValue("fm2z2");
		fm3z1 = pd->findParameterValue("fm3z1");
		fm1z3 = pd->findParameterValue("fm1z3");
		fs0 = pd->findParameterValue("fs0");
		fs1 = pd->findParameterValue("fs1");
		fs2 = pd->findParameterValue("fs2");
		fs3 = pd->findParameterValue("fs3");
		fs4 = pd->findParameterValue("fs4");
		fsz1 = pd->findParameterValue("fsz1");
		fsz2 = pd->findParameterValue("fsz2");
		fsz3 = pd->findParameterValue("fsz3");
		fsz4 = pd->findParameterValue("fsz4");
		fs1z1 = pd->findParameterValue("fs1z1");
		fs1z2 = pd->findParameterValue("fs1z2");
		fs2z1 = pd->findParameterValue("fs2z1");
		fs2z2 = pd->findParameterValue("fs2z2");
		fs3z1 = pd->findParameterValue("fs3z1");
		fs1z3 = pd->findParameterValue("fs1z3");
		p0 = pd->findParameterValue("p0");
		p1 = pd->findParameterValue("p1");
		p2 = pd->findParameterValue("p2");
		p3 = pd->findParameterValue("p3");
		p4 = pd->findParameterValue("p4");
		pz1 = pd->findParameterValue("pz1");
		pz2 = pd->findParameterValue("pz2");
		pz3 = pd->findParameterValue("pz3");
		pz4 = pd->findParameterValue("pz4");
		p1z1 = pd->findParameterValue("p1z1");
		p1z2 = pd->findParameterValue("p1z2");
		p2z1 = pd->findParameterValue("p2z1");
		p2z2 = pd->findParameterValue("p2z2");
		p3z1 = pd->findParameterValue("p3z1");
		p1z3 = pd->findParameterValue("p1z3");

		if (cm0 == -1) cm0 = 0;
		if (cm1 == -1) cm1 = 0;
		if (cm2 == -1) cm2 = 0;
		if (cm3 == -1) cm3 = 0;
		if (cm4 == -1) cm4 = 0;
		if (cmz1 == -1) cmz1 = 0;
		if (cmz2 == -1) cmz2 = 0;
		if (cmz3 == -1) cmz3 = 0;
		if (cmz4 == -1) cmz4 = 0;
		if (cm1z1 == -1) cm1z1 = 0;
		if (cm1z2 == -1) cm1z2 = 0;
		if (cm2z1 == -1) cm2z1 = 0;
		if (cm2z2 == -1) cm2z2 = 0;
		if (cm1z3 == -1) cm1z3 = 0;
		if (cm3z1 == -1) cm3z1 = 0;
		if (cs0 == -1) cs0 = 0;
		if (cs1 == -1) cs1 = 0;
		if (cs2 == -1) cs2 = 0;
		if (cs3 == -1) cs3 = 0;
		if (cs4 == -1) cs4 = 0;
		if (csz1 == -1) csz1 = 0;
		if (csz2 == -1) csz2 = 0;
		if (csz3 == -1) csz3 = 0;
		if (csz4 == -1) csz4 = 0;
		if (cs1z1 == -1) cs1z1 = 0;
		if (cs1z2 == -1) cs1z2 = 0;
		if (cs2z1 == -1) cs2z1 = 0;
		if (cs2z2 == -1) cs2z2 = 0;
		if (cs1z3 == -1) cs1z3 = 0;
		if (cs3z1 == -1) cs3z1 = 0;
		if (fm0 == -1) fm0 = 0;
		if (fm1 == -1) fm1 = 0;
		if (fm2 == -1) fm2 = 0;
		if (fm3 == -1) fm3 = 0;
		if (fm4 == -1) fm4 = 0;
		if (fmz1 == -1) fmz1 = 0;
		if (fmz2 == -1) fmz2 = 0;
		if (fmz3 == -1) fmz3 = 0;
		if (fmz4 == -1) fmz4 = 0;
		if (fm1z1 == -1) fm1z1 = 0;
		if (fm1z2 == -1) fm1z2 = 0;
		if (fm2z1 == -1) fm2z1 = 0;
		if (fm2z2 == -1) fm2z2 = 0;
		if (fm1z3 == -1) fm1z3 = 0;
		if (fm3z1 == -1) fm3z1 = 0;
		if (fs0 == -1) fs0 = 0;
		if (fs1 == -1) fs1 = 0;
		if (fs2 == -1) fs2 = 0;
		if (fs3 == -1) fs3 = 0;
		if (fs4 == -1) fs4 = 0;
		if (fsz1 == -1) fsz1 = 0;
		if (fsz2 == -1) fsz2 = 0;
		if (fsz3 == -1) fsz3 = 0;
		if (fsz4 == -1) fsz4 = 0;
		if (fs1z1 == -1) fs1z1 = 0;
		if (fs1z2 == -1) fs1z2 = 0;
		if (fs2z1 == -1) fs2z1 = 0;
		if (fs2z2 == -1) fs2z2 = 0;
		if (fs1z3 == -1) fs1z3 = 0;
		if (fs3z1 == -1) fs3z1 = 0;
		if (p0 == -1) p0 = 0;
		if (p1 == -1) p1 = 0;
		if (p2 == -1) p2 = 0;
		if (p3 == -1) p3 = 0;
		if (p4 == -1) p4 = 0;
		if (pz1 == -1) pz1 = 0;
		if (pz2 == -1) pz2 = 0;
		if (pz3 == -1) pz3 = 0;
		if (pz4 == -1) pz4 = 0;
		if (p1z1 == -1) p1z1 = 0;
		if (p1z2 == -1) p1z2 = 0;
		if (p2z1 == -1) p2z1 = 0;
		if (p2z2 == -1) p2z2 = 0;
		if (p1z3 == -1) p1z3 = 0;
		if (p3z1 == -1) p3z1 = 0;

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
	flabel = sd->findParameterValue("flabel");
	datadir = sd->findParameterValue("datadir");
	halofile = sd->findParameterValue("halofile");
	rnn_halofile = sd->findParameterValue("rnn_halofile");
	simulationfile = sd->findParameterValue("simulationfile");
	rnnfile = sd->findParameterValue("rnnfile");
	out_path = sd->findParameterValue("out_path");
	denspdffile = sd->findParameterValue("denspdffile");
	lbcgfile = sd->findParameterValue("lbcgfile");
	colortrdir = sd->findParameterValue("colortrdir");
#ifdef SHAM_TEST
	sham_file = sd->findParameterValue("SHAM_file");
#elif BCC
	PSTR = sd->findParameterValue("pstr");
	ZSTR = sd->findParameterValue("zstr");
#endif
#ifdef FITS_GALAXIES
	ffile = sd->findParameterValue("fitsfile");
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
