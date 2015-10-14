#include "utils.h"

void read_data(string filename, data_point data[], int npoints){
  ifstream infile(filename.c_str());
  if(infile.fail()){
    cerr<<"Cannot open "<<filename<<endl;
    string com = "ls -l "+filename;
    system(com.c_str());
    exit(0);
  }
  else
    for(int i=0; i<npoints; i++)
      infile>>data[i].x>>data[i].y>>data[i].err;
}

void read_cdata(string filename, vector < vector <float> > &data, int npoints){
  ifstream infile(filename.c_str());
  if(infile.fail()){
    cerr<<"Cannot open "<<filename<<endl;
    string com = "ls -l "+filename;
    system(com.c_str());
    exit(0);
  }
  else
    // cout<<"reading"<<filename<<endl;
    for(int i=0; i<npoints; i++)
      for(int j=0; j<npoints; j++)
	infile>>data[i][j];
}

/* ----------------   Splines and spline derivatives -------------------  */

void	spline(float x[], float y[] ,int n, float yp1, float ypn, float y2[])
{
int	i,k;
float	p,qn,sig,un,*u;

  if (x[0] >= x[n-1]) {
    printf("x must be monotonically increasing in spline.\n");
    exit(1);
  }
  u = (float*)malloc(n*sizeof(float));
  if (u==NULL) {perror("malloc");exit(1);}
  if (yp1 > 0.99e30)
    y2[0]=u[0]=0.0;
  else {
    y2[0] = -0.5;
    u[0]  = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
  }
  for (i=1;i<n-1;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  if (ypn > 0.99e30)
    qn=un=0.0;
  else {
    qn=0.5;
    un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  }
  y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
  for (k=n-2;k>=0;k--)
    y2[k]=y2[k]*y2[k+1]+u[k];
  free(u);
}



void	splint(float xa[], float ya[], float y2a[], int n, float x, float *y)
/// Returns the spline interpolation to the function 
{
int	klo,khi,k;
float	h,b,a;

  if (x<xa[0] || x>xa[n-1]) {
    printf("x=%e, out of range xa[0]=%e, xa[%d]=%e\n",x,xa[0],n-1,xa[n-1]);
    exit(1);
  }
  klo=0;
  khi=n-1;
  while (khi-klo > 1) {
  k=(khi+klo) >> 1;
  if (xa[k] > x) khi=k;
  else klo=k;
  }
  h=xa[khi]-xa[klo];
  if (h == 0.0) { printf("Bad xa input to routine splint\n"); exit(1); }
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}



void	spldint(float xa[], float ya[], float y2a[], int n, float x, float *y)
/// Returns the derivative of the function! 
{
int	klo,khi,k;
float	h,b,a,c;

  if (x<xa[0] || x>xa[n-1]) {
    printf("x=%e, out of range xa[0]=%e, xa[%d]=%e\n",x,xa[0],n-1,xa[n-1]);
    exit(1);
  }
  klo=0;
  khi=n-1;
  while (khi-klo > 1) {
  k=(khi+klo) >> 1;
  if (xa[k] > x) khi=k;
  else klo=k;
  }
  h=xa[khi]-xa[klo];
  if (h == 0.0) { printf("Bad xa input to routine spldint\n"); exit(1); }
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
 c=(ya[khi]-ya[klo])/h;
  *y=c-(3*a*a-1)/6*h*y2a[klo]+(3*b*b-1)/6*h*y2a[khi];
}


/* ---------------- Mathematical functions. ------------------- */

double	gammafn2(float a, float xmin)
/* The incomplete Gamma function.  Do by integration. */
{
int	k,wt,N=5000;
double	xx,hh,xmax,sum;
  xmax = 10;
  hh   = (xmax-xmin)/N;
  if (xmin>0 && xmin<xmax) {
    sum  = pow(xmin,a-1)*exp(-xmin)+0;
    for (wt=4, k=1; k<N; k++) {
      xx   = xmin + k*hh;
      sum += pow(xx,(double) (a-1))*exp(-xx)*wt;
      wt   = 8/wt;
    }
    sum *= hh/3;
  }
  else
    sum = 0;
  return(sum);
}



double	gammln(float xx)
/* Returns ln Gamma(x) if x>0 */
{
int	j;
double	x,y,tmp,ser;
double	cof[6]={76.18009172947146,-86.50532032941677,
        24.01409824083091,-1.231739572450155,
        0.1208650973866179e-2,-0.5395239384953e-5};
  if (xx<0) {
    printf("x=%e < 0 not allowed in gammln.\n",xx);
    exit(1);
  }
  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) {
    y   += 1.0;
    ser += cof[j]/y;
  }
  return( -tmp+log(2.5066282746310005*ser/x) );
}


double	gammafn1(float z)
/* The gamma function. */
{
float	y,tmp;
  if (z >= 0.1)
    tmp = exp( gammln(z) );
  else {
    y   = 1-z;
    tmp = M_PI*y/sin(M_PI*y)/exp(gammln(2.0-z)); 
  }
  return(tmp);
}


