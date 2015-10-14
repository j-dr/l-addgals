///////////////////////////////////////////////////////////////////////
//
// computes auto-correlation function
// using a PM grid 
// assumes a cubical box
// covar3 0.08 18 12 375.0 0 375.0 1 xyzxyz.dat a 1 1 xyzxyz.dat a 1 1//////////////////////////////////////////////////////////////////////
#include "cf.h"

int cf(string infn, string outfn){
  //commented out writing of cfdata.
  const char * datafilename = infn.c_str();
  const char * outfilename = outfn.c_str();

  //initialize stuff
  for (int i=0; i<NC; i++) {
    for (int j=0; j< NC; j++) {
      for (int k=0; k<NC; k++) {
        mesh[i][j][k] =NULL;
      }
    }
  }
  int ir;
  for (ir=1; ir<= RBINS; ir++) 
    pairs[ir] = 0.0;
  max_scale = rmax;
  dlogr = (log10(rmax)-log10(rmin))/RBINS;
  int ndata = load_data(datafilename);
  cout<<ndata<<" "<<sim.Boxsize()<<endl;
  compute_pairs(pairs);
  double ndens = ndata/sim.CubeVolume(); //number density

  //output the results 
  ofstream outfile;
  outfile.open(outfilename);

  for (ir=1; ir<=RBINS; ir++) {
    double logr = log10(rmin) + (ir-0.5)*dlogr; //middle of bin, for xi
    double r = pow(10.0, logr);
    double dr = pow(10.0, log10(rmin)+ir*dlogr)-pow(10.0, log10(rmin)+(ir-1)*dlogr);
    double Np = 4.0*Pi*ndata*ndens*r*r*dr;

    double xi = pairs[ir]/Np - 1.0;
    double err =  xi*sqrt(1./pairs[ir]+2./ndata);
     cout<<r<<" "<<xi<<" "<<err<<"\n";
    if(xi > 0){
      outfile<<r<<" "<<xi<<" "<<err<<"\n";
    }
  }
  outfile.close();

  for (int i=0; i<NC; i++) {
    for (int j=0; j< NC; j++) {
      for (int k=0; k<NC; k++) {
        mesh[i][j][k] = NULL;
        delete mesh[i][j][k];
      }
    }
  }
  return 1;
}

int load_data(const char* filename) {
  data* last[NC][NC][NC];
  ifstream infile;
  infile.open(filename);
  if (infile.fail()) {
    cerr<<"cannot open file "<<filename<<"\n";
    exit(1);
  }
#ifndef QUIET
  cout<<"#reading "<<filename<<"...\n"; cout.flush();
#endif
  char buff;
  char lastbuff;
  int ndata=0;
  int ntot=0;
  //infile>>Lbox;
  //cellsize = Lbox/NC;
#ifndef QUIET
  cout<<"cellsize: "<<cellsize<<" max scale: "<<max_scale<<"\n";
#endif
  //  int halo_id=0;
  while(infile>>buff) {
    lastbuff=buff;
    infile.putback(buff);
    ntot++;
    data* p = new data;   
    //cout<<"%%%%%%%%%"<<p<<endl;
    int keep;
    keep = p->get_points_data(&infile);
    if (keep) {
      p->attach(last);  
      ndata++;
    }
    else
      delete p;
  } //not eof  
#ifndef QUIET 
  cout<<"#last halo: "<<halo_id<<"\n";
#endif
  infile.close();
  //cout<<ndata<<" points read\n";
#ifndef QUIET 
  cout<<ntot<<" halos read, "<<ndata<<" halos/galaxies kept ("
      <<((float)ndata)/((float)ntot)<<")\n";
#endif
  //(ndata>maxnum) cullfrac = cullfrac*maxnum/ndata;   
  return ndata;
}




void data::attach(data* last[NC][NC][NC]) {
  next = NULL;
  int i = (int)(pos.x/cellsize);
  if (i==NC) i=NC-1;
  int j = (int)(pos.y/cellsize);
  if (j ==NC) j=NC-1;
  int k = (int)(pos.z/cellsize);
  if (k ==NC) k=NC-1;
  if (i >= NC || j >= NC || k >= NC) {
    cerr<<"LBOX: "<<sim.Boxsize()<<" "<<pos.x<<" "<<pos.y<<" "<<pos.z<<"\n";
    cerr<<"point out of bounds: check box size\n";
    exit(1);
  }
  // p->pos.print(); cout<<" : "<<i<<" "<<j<<" "<<k<<"\n";
  if (mesh[i][j][k] == NULL) {
    //cout<<"new point\n";
    mesh[i][j][k] = this;
    last[i][j][k] = this;
  }
  else {
    //cout<<"adding to "; last[i][j][k]->pos.print(); cout<<"\n";
    last[i][j][k]->next = this;
    last[i][j][k] = this;
  }
}


void compute_pairs(double* pairs) {
  //loop over every point in box, looking for pairs
  int n=0;
#ifndef QUIET
  cout<<"computing pairs...\n"; cout.flush();
#endif
  for (int i1=0; i1 < NC; i1++) {
    for (int j1=0; j1 < NC; j1++) {
      for (int k1=0; k1 < NC; k1++) {
	data* p1 = mesh[i1][j1][k1];
	while (p1) { //first member of pair
	  n++;
	  //find region of grid to search for neighbors
	  
	  int il, jl, kl, ih, jh, kh; //regular indices
	  int ilp, jlp, klp, ihp, jhp, khp; //periodic indices
	  //set up limits for periodic boundary conditions
	  int pbc_xl = setup_lower_index(p1->pos.x, il, ilp); 
	  int pbc_yl = setup_lower_index(p1->pos.y, jl, jlp);
	  int pbc_zl = setup_lower_index(p1->pos.z, kl, klp);
	  int pbc_xh = setup_higher_index(p1->pos.x, ih, ihp); 
	  int pbc_yh = setup_higher_index(p1->pos.y, jh, jhp);
	  int pbc_zh = setup_higher_index(p1->pos.z, kh, khp);
	  //search grid
	  vec3 wrap; //converts coordinates properly using periodic bc
	  //set up starting indeces for x loop
	  int xloop=1; int i2; looptypetype xlooptype;
	  setup_firstloop(pbc_xl, il, ilp, i2, xlooptype, wrap.x);
	  while (xloop) {
	    //set up starting indeces for y loop	
	    int yloop=1; int j2; looptypetype ylooptype;
	    setup_firstloop(pbc_yl, jl, jlp, j2, ylooptype, wrap.y);
	    while (yloop) {
	      //set up starting indeces for z loop
	      int zloop=1; int k2; looptypetype zlooptype;
	      setup_firstloop(pbc_zl, kl, klp, k2, zlooptype, wrap.z);
	      while (zloop) {
		vec3 cell_center((i2+0.5)*cellsize, (j2+0.5)*cellsize, 
				   (k2+0.5)*cellsize);
		vec3 d12 = cell_center+wrap-p1->pos;
		double min_distance = d12.magnitude()-cellsize*sqrt(3.0)/2.0;
		if (min_distance < max_scale) {
		  data* p2 = mesh[i2][j2][k2];
		  while (p2) { //step through linked list of neighbors
		    vec3 wrapped_pos2 = p2->pos+wrap;
		    vec3 r = p1->pos-wrapped_pos2;
		    //this is the magnitude of the vector r, not some galaxy!
		    double mag_r = log10(r.magnitude());
		    cout<<mag_r<<" "<<log10(wrapped_pos2.magnitude())<<endl;
		    if (mag_r > log10(rmin) && mag_r <=log10(rmax)) {
		      int rbin = (int)((mag_r - log10(rmin))/dlogr)+1;
		      //*********************************************
		      //*********************************************
		      if(paircount==ALL){
			if (rbin >=0 && rbin <= RBINS) 
			  pairs[rbin] += 1.0; //each pair is counted twice
		      }
		      else if(paircount==NOSAMEHALO){
			if(p1->pos.id!=p2->pos.id)			
			  if (rbin >=0 && rbin <= RBINS) 
			    pairs[rbin] += 1.0; //each pair is counted twice
		      }
		      else if(paircount==SAMEHALO){
			if(p1->pos.id==p2->pos.id)			
			  if (rbin >=0 && rbin <= RBINS) 
			    pairs[rbin] += 1.0; //each pair is counted twice
		      }
		    }
		    p2 = p2->next;
		  } //while p2
		} //mindistance < max_scale
		k2++; //go on to next cell
		zloop=setup_nextloop(pbc_zh, kh, khp, k2, zlooptype,wrap.z);
	      } //k2 (zloop)
	      j2++;
	      yloop=setup_nextloop(pbc_yh, jh, jhp, j2, ylooptype, wrap.y);
	    } //j2 (yloop)
	    i2++;
	    xloop=setup_nextloop(pbc_xh, ih, ihp, i2, xlooptype, wrap.x);
	  } //i2 (xloop)
	  p1 = p1->next;
	} //while p1
      } //k1
    } //j1
  } //i1
}

int setup_lower_index(double x, int &il, int &ilp) { 
  double cell = (x-max_scale)/cellsize;
  int pbc;
  if (cell < 0.0) { //deal with periodic boundary conditions
    ilp = (int)(NC+cell);
    il = 0;
    pbc =1;
  }
  else {
    il = (int)cell;
    ilp = -1; 
    pbc=0;
  }
  return pbc;
}

int setup_higher_index(double x, int &ih, int &ihp) { 
  double cell = (x+max_scale)/cellsize;
  int pbc;
  if (cell > NC-1) {
    ihp = (int)(cell-NC);
    ih = NC-1;
    pbc =1;
  }
  else {
    ih = (int)cell;
    ihp = -1;
    pbc=0;
  }
  return pbc;
}

void setup_firstloop(int pbc, int il, int ilp, 
		     int &i2, looptypetype &looptype, double &wrap) {
  if (pbc) {
    i2 = ilp;
    wrap = -sim.Boxsize();
    looptype = PERIODIC_L;
  }
  else {
    i2 = il;
    wrap = 0.0;
    looptype = NORMAL;
  }
}


int setup_nextloop(int pbc, int ih, int ihp, int &i2,
		    looptypetype &looptype, double &wrap) {
  int loop=1;
  if (looptype == PERIODIC_L && i2 > NC-1) { 
    i2 = 0;
    looptype = NORMAL;
    wrap = 0.0;
  } 
  else if (looptype == NORMAL && i2 > ih) {
    if (pbc) {
      i2 = 0;
      looptype = PERIODIC_H;
      wrap = sim.Boxsize();
    }
    else
      loop = 0; //stop loop
  }
  else if (looptype == PERIODIC_H && i2 > ihp) 
    loop =0;
  return loop;
}

