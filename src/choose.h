#ifndef choose_h
#define choose_h
#endif
#include "particle.h"
#include "galaxy.h"
#include "color.h"
#include "stl_util.h"
#include "myrand.h"

using std::cout;

//#define NODENSCORR
//#define NOMAGCORR

class CloseParticle{
public:
  CloseParticle(float value, float within):_val(value),_within(within){}
  float  comp_val() const {return _val;}
  float slop_val() const {return _within;}
  //  bool operator()(Particle p) const;
  bool operator()(Particle *p) const;
  
private:
  float _val;
  float _within;
};

inline bool CloseParticle::operator()(Particle *p) const {
  float logvalue=log10(p->Dist8());
  //return ((value>_val-_within)&&(value<_val+_within));
return ((logvalue>_val-_within)&&(logvalue<_val+_within));
}

class CloseGalaxy{
public:
  CloseGalaxy(float value, float within):_val(value),_within(within){}
  float  comp_val() const {return _val;}
  float slop_val() const {return _within;}
  bool operator()(Galaxy *g) const;
  
private:
  float _val;
  float _within;
};

inline bool CloseGalaxy::operator()(Galaxy *g) const {
  float logvalue=g->Dist8();
  return ((logvalue>_val-_within)&&(logvalue<_val+_within));
}

#ifdef SEDDLUM

class CloseSED{
public:
  CloseSED(float mag, float dens, float dens_within):
  _mag(mag),_dens(dens),_dens_within(dens_within){
    float vol_const = 1.65e6;
    //cout<<"geting n "<<mag<<endl;
    float this_dens = LumNumberDensityInterp(mag);
    //cout<<"Densities for Mr = -22., -18., -14., -10., -8., -4."<<endl;
    //cout<<"           "<<LumNumberDensityInterp(-22.)<<", "<<LumNumberDensityInterp(-18.)<<", "<<LumNumberDensityInterp(-14.)<<", "<<LumNumberDensityInterp(-10.)<<", "<<LumNumberDensityInterp(-8.)<<", "<<LumNumberDensityInterp(-4.)<<endl;
    //cout<<" check: this_dens = "<<this_dens<<" for magnitude "<<mag<<endl;;
    if (!((this_dens >= LumNumberDensityInterp(-27)) && (this_dens <= LumNumberDensityInterp(-1.0)))){
      cout<<"We've got a bad dens value!  Going to fail assert soon."<<endl;
      cout<<LumNumberDensityInterp(-27)<<" "<<this_dens<<" "<<LumNumberDensityInterp(-10.0)<<endl;
      cout<<LumNumberDensityInterp(Magmin)<<endl;
      cout<<"This came from a magnitude of "<<mag<<endl;
    }
    float nhere = this_dens*vol_const;
    float nddown, ndup;
    if(nhere>51){
      nddown = (nhere-50)/vol_const;
      ndup = (nhere+50)/vol_const;
    }
    else{
      nddown = (.8)/vol_const;  //just to make sure you include them all
      ndup = (101.)/vol_const;
    }
    //    cout<<mag<<" "<<nddown<<" "<<ndup<<" ";
    if (!((this_dens >= LumNumberDensityInterp(-27)) && (this_dens <= LumNumberDensityInterp(-1.0))))
      cout<<"Starting NdensMagnitude for bad dens with ndup = "<<ndup<<endl;
    _max_mag = NdensMagnitude(ndup);
    if (!((this_dens >= LumNumberDensityInterp(-27)) && (this_dens <= LumNumberDensityInterp(-1.0))))
      cout<<"Starting NdensMagnitude for bad dens with ndown = "<<nddown<<endl;
    _min_mag = NdensMagnitude(nddown);
    //cout<<_min_mag<<" "<<_max_mag<<" ";

    //    if(_min_mag<Magmin) _min_mag = Magmin-0.01;

    float min_tol = 0.05;
    if(mag > -19.) min_tol = 0.5;
    if(mag > -13.5) min_tol = 1.0;
    if(mag > -11.) min_tol = 2.0;
    if(mag > -10.) min_tol = 3.0;
    if(mag > -9.) min_tol = 4.0;
    if(mag > -8.) min_tol = 5.0;
    if(mag > -7.) min_tol = 6.0;

    if(_min_mag>mag-min_tol){
      //     cout<<mag<<" min "<<_min_mag<<" "<<nddown<<endl;
      _min_mag = mag-min_tol;
    }
    if(_max_mag<mag+min_tol) {
      //cout<<mag<<
      _max_mag = mag+min_tol;
    }

    //hack added by mbusha to deal with dim galaxies
    //b/c of the mag distribution, we'll at least select from galaxies with -19<Mr<0
    //float mag_cut = -21.0;
    float mag_cut = -20.0;
    if(_min_mag > mag_cut) _min_mag = mag_cut;

    //hack by mbusha to deal w/ fact that dr3 catalog only has galaxies with mr>-22.34
    //cout<<"_min_mag = "<<_min_mag<<" _max_mag = "<<_max_mag<<" mag = "<<mag<<endl;
    float brightest_sed = -22.0;
    if(_max_mag < brightest_sed) _max_mag = brightest_sed;
    //cout<<"_min_mag = "<<_min_mag<<" _max_mag = "<<_max_mag<<" mag = "<<mag<<endl;
	

    //    cout<<min_mag<<" "<<mag<<" "<<max_mag<<endl;
    //modified 06-06-04
  }
  bool operator()(GalSED g) const;
  
private:
  float _mag;
  //  float _mag_within;
  float _max_mag;
  float _min_mag;

  float _dens;
  float _dens_within;
};

inline bool CloseSED::operator()(GalSED g) const {
  float mag=g.MR();
  float logdens=log10(g.Dens());
  //  if((_dens==0)||_mag<-21.9))
  if((_dens==0))
    return ((mag>_min_mag)&&(mag<_max_mag));
    // changed 03-02-04
    // return ((mag>_mag-_mag_within)&&(mag<_mag+_mag_within)&&
    //	    (dens==_dens));
  else{
#ifndef COLORS_FROM_RELATIVE_DENSITY
    logdens=log10(g.Dens());
    float _logdens=log10(_dens);
#else
    logdens=g.Dens();
    float _logdens=_dens;
#endif
    return((mag>_min_mag)&&(mag<_max_mag)&&
	   (logdens>(_logdens-_dens_within))&&(logdens<(_logdens+_dens_within)));
  }
}

class CloseSEDMag{
public:
  CloseSEDMag(float value, float within):_val(value),_within(within){}
  float  comp_val() const {return _val;}
  float slop_val() const {return _within;}
  bool operator()(GalSED g) const;
  
private:
  float _val;
  float _within;
};

inline bool CloseSEDMag::operator()(GalSED g) const {
  float value=g.MR();
  return ((value>_val-_within)&&(value<_val+_within));
}

class CloseSEDDens{
public:
  CloseSEDDens(float value, float within):_val(value),_within(within){}
  float  comp_val() const {return _val;}
  float slop_val() const {return _within;}
  bool operator()(GalSED g) const;
  
private:
  float _val;
  float _within;
};

inline bool CloseSEDDens::operator()(GalSED g) const {
  float value=g.Dens();
  if(_val==0)
    return    (value==_val);
  else return
	 ((log10(value)>log10(_val)-_within)&&(log10(value)<log10(_val)+_within));

}

// choose a color whose mag and dens are close to what you are looking for.
inline int ChooseSED(vector <GalSED> v, 
		     double mag, double dens, double ThisZ, int ThisBCG){

  int mag_min = -22.2;
  //if (mag < mag_min)
  cout<<"Looking at galaxy with mr = "<<mag<<" and dens = "<<dens<<endl;

  float dens_within = 0.12;
  //  cout<<"dens="<<dens<<" "<<"mag"<<mag<<endl;
  if(dens >0){
    double tmp = log10(dens);
    dens_within = 0.035;
    
    // from choose.pro
  // histogauss, alog10(d10), a, x, y
  // getbins, x, y, bin, t
    
    if(tmp<    -0.784051) dens_within =      0.369258;
    else if(tmp<    -0.765924) dens_within =      0.376225;
    else if(tmp<    -0.747797) dens_within =      0.371551;
    else if(tmp<    -0.729670) dens_within =      0.374431;
    else if(tmp<    -0.711542) dens_within =      0.352473;
    else if(tmp<    -0.693415) dens_within =      0.356601;
    else if(tmp<    -0.675288) dens_within =      0.333629;
    else if(tmp<    -0.657161) dens_within =      0.332235;
    else if(tmp<    -0.639034) dens_within =      0.309356;
    else if(tmp<    -0.620907) dens_within =      0.294840;
    else if(tmp<    -0.602779) dens_within =      0.285336;
    else if(tmp<    -0.584652) dens_within =      0.266087;
    else if(tmp<    -0.566525) dens_within =      0.254914;
    else if(tmp<    -0.548398) dens_within =      0.240171;
    else if(tmp<    -0.530271) dens_within =      0.226590;
    else if(tmp<    -0.512143) dens_within =      0.213583;
    else if(tmp<    -0.494016) dens_within =      0.200794;
    else if(tmp<    -0.475889) dens_within =      0.195970;
    else if(tmp<    -0.457762) dens_within =      0.177507;
    else if(tmp<    -0.439635) dens_within =      0.167844;
    else if(tmp<    -0.421507) dens_within =      0.153428;
    else if(tmp<    -0.403380) dens_within =      0.156369;
    else if(tmp<    -0.385253) dens_within =      0.146977;
    else if(tmp<    -0.367126) dens_within =      0.129867;
    else if(tmp<    -0.348999) dens_within =      0.125221;
    else if(tmp<    -0.330872) dens_within =      0.126469;
    else if(tmp<    -0.312744) dens_within =      0.113973;
    else if(tmp<    -0.294617) dens_within =      0.111962;
    else if(tmp<    -0.276490) dens_within =      0.113973;
    else if(tmp<    -0.258363) dens_within =      0.100089;
    else if(tmp<    -0.240236) dens_within =     0.0888588;
    else if(tmp<    -0.222108) dens_within =     0.0879961;
    else if(tmp<    -0.203981) dens_within =     0.0799730;
    else if(tmp<    -0.185854) dens_within =     0.0740894;
    else if(tmp<    -0.167727) dens_within =     0.0799730;
    else if(tmp<    -0.149600) dens_within =     0.0684045;
    else if(tmp<    -0.131473) dens_within =     0.0623285;
    else if(tmp<    -0.113345) dens_within =     0.0617973;
    else if(tmp<   -0.0952182) dens_within =     0.0589503;
    else if(tmp<   -0.0770910) dens_within =     0.0592715;
    else if(tmp<   -0.0589638) dens_within =     0.0567954;
    else if(tmp<   -0.0408366) dens_within =     0.0562084;
    else if(tmp<   -0.0227094) dens_within =     0.0459563;
    else if(tmp<  -0.00458223) dens_within =     0.0418320;
    else if(tmp<    0.0135449) dens_within =     0.0310752;
    else if(tmp<    0.0316721) dens_within =     0.0286219;
    else if(tmp<    0.0497993) dens_within =     0.0327600;
    else if(tmp<    0.0679265) dens_within =     0.0356601;
    else if(tmp<    0.0860537) dens_within =     0.0281770;
    else if(tmp<     0.104181) dens_within =     0.0297981;
    else if(tmp<     0.122308) dens_within =     0.0308111;
    else if(tmp<     0.140435) dens_within =     0.0261450;
    else if(tmp<     0.158562) dens_within =     0.0267232;
    else if(tmp<     0.176690) dens_within =     0.0271230;
    else if(tmp<     0.194817) dens_within =     0.0255913;
    else if(tmp<     0.212944) dens_within =     0.0243864;
    else if(tmp<     0.231071) dens_within =     0.0166815;
    if(tmp>     0.249198) dens_within =     0.0275350;
    if(tmp>     0.267326) dens_within =     0.0177718;
    if(tmp>     0.285453) dens_within =     0.0257733;
    if(tmp>     0.303580) dens_within =     0.0250606;
    if(tmp>     0.321707) dens_within =     0.0243864;
    if(tmp>     0.339834) dens_within =     0.0245515;
    if(tmp>     0.357961) dens_within =     0.0250606;
    if(tmp>     0.376089) dens_within =     0.0308111;
    if(tmp>     0.394216) dens_within =     0.0269216;
    if(tmp>     0.412343) dens_within =     0.0279597;
    if(tmp>     0.430470) dens_within =     0.0175425;
    if(tmp>     0.448597) dens_within =     0.0397914;
    if(tmp>     0.466725) dens_within =     0.0330587;
    if(tmp>     0.484852) dens_within =     0.0399865;
    if(tmp>     0.502979) dens_within =     0.0302962;
    if(tmp>     0.521106) dens_within =     0.0431600;
    if(tmp>     0.539233) dens_within =     0.0445751;
    if(tmp>     0.557360) dens_within =     0.0479838;
    if(tmp>     0.575488) dens_within =     0.0533153;
    if(tmp>     0.593615) dens_within =     0.0639784;
    if(tmp>     0.611742) dens_within =     0.0637907;
    if(tmp>     0.629869) dens_within =     0.0804462;
    if(tmp>     0.647996) dens_within =     0.0847066;
    if(tmp>     0.666124) dens_within =      0.102930;
    if(tmp>     0.68) dens_within =      0.12;
    if(tmp>     0.70) dens_within =      0.13;
    if(tmp>     0.75) dens_within =      0.14;
    if(tmp>     0.80) dens_within =      0.15;
    if(tmp>     0.85) dens_within =      0.16;
    if(tmp>     0.90) dens_within =      0.17;
    if(tmp>     0.95) dens_within =      0.18;
    if(tmp>     1.00) dens_within =      0.19;
    //if dens is zero then it looks for other zero elements
    // dens_within = 0.5*dens_within;
  }
#ifdef COLORS_FROM_RELATIVE_DENSITY
#ifdef RED_FRACTION
  dens_within = 0.0025;
#else
  dens_within = 0.0005; //typically 10-20 choices
#endif
  //dens_within = 0.0002; //typically 2-10 choices
#endif
#ifndef NODENSCORR
#ifndef NOMAGCORR  
  // cout<<"using"<<dens<<endl;
  //cout<<"Making constraing with mag = "<<mag<<", dens = "<<dens<<", dens_within = "<<dens_within<<endl;
  //if (mag < mag_min)
  //cout<<"Finding SED, dens_within = "<<dens_within<<endl;
  CloseSED cc(mag, dens, dens_within);
  //if (mag < mag_min)
  //cout<<"Some SEDs located"<<endl;
#endif
#endif
#ifdef NODENSCORR
  CloseSEDMag cc(mag, mag_within);
#endif
#ifdef NOMAGCORR
  CloseSEDDens cc(dens, dens_within);
#endif

  vector <GalSED> tmp_v;
  tmp_v.reserve(50);
  // copy all the galaxies which are close enough
  //cout<<"copy_if"<<endl;
  copy_if(v.begin(),v.end(),
          back_inserter(tmp_v), cc);
  //if (mag < mag_min)
  //cout<<"Copied "<<tmp_v.size()<<" SEDs"<<endl;
  //cout<<"copied "<<tmp_v.size()<<" sed's."<<endl;
  // ------> USE THIS LINE TO SEE HOW MANY CANDIDATE GALAXIES DENS_WITHIN CHOOSES <------------ \\
  //cout<<tmp_v.size()<<" "<<mag<<" "<<dens<<" "<<dens_within<<endl;
  //if(tmp_v.size()==0){ cout<<tmp_v.size()<<" tmp_vsize:" <<mag<<" "<<dens<<" ";}// system("date");};
  //cout<<tmp_v.size()<<" ";
    //just in case you didn't find any -- this should only happen rarely
  while(tmp_v.size()<1){
    // copy all the particles which are close enough
    //    mag_within = mag_within+0.1;
    if(dens > 0)
      dens_within = dens_within+0.04;
    else{
      cout<<mag<<" "<<dens<<" "<<dens_within<<endl;  
      cout<<"Infinite loop... exiting"<<endl;
    }
    //if (mag < mag_min)
    //cout<<"expanded dens_within = "<<dens_within<<", dens = "<<dens<<", mag = "<<mag<<endl;
    CloseSED cc2(mag, dens, dens_within);
    copy_if(v.begin(),v.end(),
            back_inserter(tmp_v), cc2);
    //cout<<"   expand to "<<tmp_v.size()<<" "<<dens_within<<endl;
    //cout<<tmp_v.size()<<" "<<-1<<" "<<dens<<" "<<dens_within<<endl;
    //   cout<<tmp_v.size()<<" "; 
    //system("date");
    //system(date);
  }
  //if (mag < mag_min)
  //cout<<"Choosing from "<<tmp_v.size()<<" SEDs"<<endl;
  int ri = randint(tmp_v.size())-1;
  //if (mag < mag_min)
  //cout<<"Chose SED "<<ri<<endl;
  //cout<<tmp_v.size()<<" "<<mag<<" "<<dens<<" "<<tmp_v[ri].MR()<<" "<<tmp_v[ri].Dens()<<" "<<endl;
#ifdef RED_FRACTION
  //choose our SED to match "red fraction observations"
  cout<<"doing red fraction for particle with z = "<<ThisZ<<" and "<<tmp_v.size()<<" candidates."<<endl;
  int nred = 0;
  for(int i=0;i<tmp_v.size();i++)
    nred += v[tmp_v[i].Id()].Red();
  int nblue = tmp_v.size() - nred;
  float local_red_fraction = ((float) nred) / ((float) tmp_v.size());
  cout<<"  nred = "<<nred<<", mblue = "<<nblue<<endl;
  cout<<"  ncandidates = "<<tmp_v.size()<<endl;
  if (nred == 0 || nred == tmp_v.size()){
    return tmp_v[ri].Id();
  }
  /*
  //float redfrac_start_evolve = 0.3;
  float redfrac_start_evolve = 0.2;
  float redfraction03 = 0.66;
  float redfraction08 = 0.14;
  //float redfraction08 = 0.12;
  float slope = (redfraction03 - redfraction08)/(redfrac_start_evolve-0.8);
  float intercept = redfraction03 - slope*redfrac_start_evolve;
  if (ThisZ > redfrac_start_evolve && ThisBCG == 0){
  */
  float slope = (REDFRACTION1 - REDFRACTION2)/(Z_REDFRACTION1-Z_REDFRACTION2);
  float intercept = REDFRACTION1 - slope*Z_REDFRACTION1;
  if (ThisZ > Z_REDFRACTION1 && ThisBCG == 0){
    cout<<"  Resetting fraction."<<endl;
    float global_fraction = slope*ThisZ + intercept;
    if (global_fraction < REDFRACTION2)
      global_fraction = REDFRACTION2;
    cout<<"   global_fraction = "<<global_fraction<<endl;
    cout<<"   local_fraction = "<<local_red_fraction<<endl;

    /*
    float ran_red = drand48();
    float val = 0.;
    int ind = -1;
    if (ran_red < global_fraction){ //select a red galaxy
      float ran = drand48()*nred;
      for(ind=0,val=0;val<tmp_v.size();val++){
        if(v[tmp_v[ind].Id()].Red() == 1)
	  val++;
	if (val >= ran) break;
      }
    } else {
      float ran = drand48()*nblue;
      for(ind=0,val=0;val<tmp_v.size();val++){
        if(v[tmp_v[ind].Id()].Red() == 0)
          val++;
        if (val >= ran) break;
      }
    }
    */

    float fac = global_fraction/REDFRACTION1;
    float target_red_fraction = local_red_fraction*fac;
    cout<<"   fac = "<<fac<<endl;
    cout<<"   target_red_fraction = "<<target_red_fraction<<endl;
    float newfac = target_red_fraction*((float) tmp_v.size() - nred) / (nred*(1.0 - target_red_fraction));
    float new_size = tmp_v.size() - nred + nred*newfac;
    cout<<"  newfac = "<<newfac<<endl;
    cout<<"  new_size = "<<new_size<<endl;
    float ran = drand48()*new_size;
    float val = 0;
    int ind = 0;
    //cout<<"   Looking for random number "<<ran<<endl;
    for(ind=0;ind<tmp_v.size();ind++){
      if (v[tmp_v[ind].Id()].Red() == 1)
	{
	  val += newfac;
	}
      else
	{
	  val++;
	}
      if (val > ran)
	break;
    }

    //cout<<"Selected: "<<ind<<", val = "<<val<<endl;
    return tmp_v[ind].Id();
  }
  else{
    return tmp_v[ri].Id();
  }
#else
  return tmp_v[ri].Id();
#endif

}

#endif

// choose a color whose mag is close to what you are looking for.
inline int ChooseBCG_SED(vector <GalSED> v, double mag){
  float dens_within = 0.12;
  //  cout<<"dens="<<dens<<" "<<"mag"<<mag<<endl;
  float dens = 0.0;

#ifndef NODENSCORR
#ifndef NOMAGCORR  
  // cout<<"using"<<dens<<endl;
  CloseSED cc(mag, dens, dens_within);
#endif
#endif
#ifdef NODENSCORR
  CloseSEDMag cc(mag, mag_within);
#endif
#ifdef NOMAGCORR
  CloseSEDDens cc(dens, dens_within);
#endif

  vector <GalSED> tmp_v;
  tmp_v.reserve(50);
  // copy all the galaxies which are close enough
  copy_if(v.begin(),v.end(),
          back_inserter(tmp_v), cc);
  //cout<<tmp_v.size()<<" "<<mag<<" "<<dens<<" "<<dens_within<<endl;
  //if(tmp_v.size()==0){ cout<<tmp_v.size()<<" tmp_vsize:" <<mag<<" "<<dens<<" ";}// system("date");};
  //cout<<tmp_v.size()<<" ";
    //just in case you didn't find any -- this should only happen rarely
  while(tmp_v.size()<1){
    // copy all the particles which are close enough
    //    mag_within = mag_within+0.1;
    if(dens > 0)
      dens_within = dens_within+0.04;
    else{
      cout<<mag<<" "<<dens<<" "<<dens_within<<endl;  
      cout<<"Infinite loop... exiting"<<endl;
    }
    CloseSED cc2(mag, dens, dens_within);
    copy_if(v.begin(),v.end(),
            back_inserter(tmp_v), cc2);
    //cout<<tmp_v.size()<<" "<<-1<<" "<<dens<<" "<<dens_within<<endl;
    //   cout<<tmp_v.size()<<" "; 
    //system("date");
    //system(date);
  }
  int ri = randint(tmp_v.size())-1;
  //cout<<tmp_v.size()<<" "<<mag<<" "<<dens<<" "<<tmp_v[ri].MR()<<" "<<tmp_v[ri].Dens()<<" "<<endl;
  return tmp_v[ri].Id();


}
