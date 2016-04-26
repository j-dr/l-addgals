#include "healpix_utils.h"
#include <iostream>
#include <cmath>
#include <climits>
#include <algorithm>
#include <functional>
#include <vector>
#include <utility>
#include <cassert>

//using namespace std;

static long ctab[0x100];
static long utab[0x100];
static const long jrll[] = { 2,2,2,2,3,3,3,3,4,4,4,4 };
static const long jpll[] = { 1,3,5,7,0,2,4,6,1,3,5,7 };
static int HEALPIX_TOOLS_INIT = 1;

struct ToPeano : public std::unary_function<long,long> {
  const long order_;
  ToPeano(long o) : order_(o) {}
  long operator() (long pix) {return nest2peano(pix, order_);}
};


long imodulo(long v1, long v2)
{ 
  return (v1>=0) ? ((v1<v2) ? v1 : (v1%v2)) : ((v1%v2)+v2);
}


long isqrt(long i)
{
  return sqrt(((double) (i)) + 0.5);
}

void tablefiller(void)
{
  assert(CHAR_BIT == 8);

  int m;
  for (m=0; m<0x100; ++m)
    {
      ctab[m] =
	(m&0x1 )       | ((m&0x2 ) << 7) | ((m&0x4 ) >> 1) | ((m&0x8 ) << 6)
	| ((m&0x10) >> 2) | ((m&0x20) << 5) | ((m&0x40) >> 3) | ((m&0x80) << 4);
      utab[m] =
	(m&0x1 )       | ((m&0x2 ) << 1) | ((m&0x4 ) << 2) | ((m&0x8 ) << 3)
	| ((m&0x10) << 4) | ((m&0x20) << 5) | ((m&0x40) << 6) | ((m&0x80) << 7);
    }
}

long nest2peano(long pix, long order_)
{
  static const unsigned long subpix[8][4] = {
    { 0, 1, 3, 2 }, { 3, 0, 2, 1 }, { 2, 3, 1, 0 }, { 1, 2, 0, 3 },
    { 0, 3, 1, 2 }, { 1, 0, 2, 3 }, { 2, 1, 3, 0 }, { 3, 2, 0, 1 } };
  static const unsigned long subpath[8][4] = {
    { 4, 0, 6, 0 }, { 7, 5, 1, 1 }, { 2, 4, 2, 6 }, { 3, 3, 7, 5 },
    { 0, 2, 4, 4 }, { 5, 1, 5, 3 }, { 6, 6, 0, 2 }, { 1, 7, 3, 7 } };
  static const unsigned long face2path[12] = {
    2, 5, 2, 5, 3, 6, 3, 6, 2, 3, 2, 3 };
  static const unsigned long face2peanoface[12] = {
    0, 5, 6, 11, 10, 1, 4, 7, 2, 3, 8, 9 };
  
  long npix_ = 1;
  npix_ = 12*(npix_ << (2*order_));
  assert(pix >= 0 && pix < npix_);
  
  long face = pix>>(2*order_);
  unsigned long path = face2path[face];
  long result = 0;
  
  long shift;
  for(shift=2*order_-2; shift>=0; shift-=2)
    {
      unsigned char spix = (pix>>shift) & 0x3;
      result <<= 2;
      result |= subpix[path][spix];
      path=subpath[path][spix];
    }

  return result + ((face2peanoface[face])<<(2*order_));
}

void ring2peanoindex(long pix, long order1_, long order2_, std::vector<long> &pidx)
{

  int nmap = 12 * ( 1 << ( 2 * order2_ ) );
  std::vector<long> hopix( 1 << 2*( order2_ - order1_ ) );

  pix = ring2nest( pix, order1_ );
  higher_nest( pix, order1_, order2_, &hopix[0] );

  /*std::cout << "Higher order pixels for nest pix " << pix << std::endl;
  for (std::vector<long>::iterator itr=hopix.begin(); itr!=hopix.end(); itr++)
    {
      std::cout << *itr << std::endl;
    }
  */
  std::transform( hopix.begin(), hopix.end(), pidx.begin(), 
  		  ToPeano(order2_) );
  //std::copy( hopix.begin(), hopix.end(), pidx.begin() );
  std::sort( pidx.begin(), pidx.end() );

  /*std::cout << "Peano indices for higher order pixels for nest pix" << pix << std::endl;
  for (std::vector<long>::iterator itr=pidx.begin(); itr!=pidx.end(); itr++)
    {
      std::cout << *itr << std::endl;
    }
  */
}

long nest2ring(long pix, long order_)
{
  long ix, iy, face_num;
  nest2xyf(pix,&ix,&iy,&face_num,order_);
  return xyf2ring(ix,iy,face_num,order_);
}

long ring2nest(long pix, long order_)
{
  long ix, iy, face_num;
  ring2xyf(pix,&ix,&iy,&face_num,order_);
  return xyf2nest(ix,iy,face_num,order_);
}

void higher_nest(long pix, long order1_, long order2_, long *hopix)
{
  assert(order2_ > order1_);
  int i;
  int base = pix >> 2 * order1_;
  int subpix = pix & ( ( 1 << ( 2 * order1_ ) ) - 1 );
  for (i = 0; i < (1 << 2 * ( order2_ - order1_ )); i++)
    {
      hopix[i] = ( ( base * ( 1 << ( 2 * order1_ ) ) + subpix ) 
	           << 2 * ( order2_ - order1_ ) ) + i;
    }
}

long lower_nest(long pix, long order1_, long order2_)
{
  assert(order1_ > order2_);
  int i;
  int base = pix >> 2 * order1_;
  int hosubpix = pix & ( ( 1 << ( 2 * order1_ ) ) - 1 );
  int losubpix = hosubpix / ( 1 << ( 2 * ( order1_ - order2_) ) );

  return base * ( 1 << ( 2 * order2_ ) ) + losubpix;
}

long xyf2nest(long ix, long iy, long face_num, long order_)
{
  if(HEALPIX_TOOLS_INIT)
    {
      tablefiller();
      HEALPIX_TOOLS_INIT = 0;
    }
  
  long pix = ((face_num)<<(2*order_)) +
    (   ((utab[ ix     &0xff]))
	| ((utab[(ix>> 8)&0xff])<<16)
	| ((utab[(ix>>16)&0xff])<<32)
	| ((utab[(ix>>24)&0xff])<<48)
	| ((utab[ iy     &0xff])<<1)
	| ((utab[(iy>> 8)&0xff])<<17)
	| ((utab[(iy>>16)&0xff])<<33)
	| ((utab[(iy>>24)&0xff])<<49) ); 
  
  long npix_ = 1;
  npix_ = 12*(npix_ << (2*order_));
  assert(pix >= 0 && pix < npix_);
  
  return pix;
}

long xyf2ring(long ix, long iy, long face_num, long order_)
{
  long npix_ = 1;
  npix_ = 12*(npix_ << (2*order_));
    
  long nside_ = 1;
  nside_ = nside_ << order_;
  
  long npface_ = 1;
  npface_ = npface_ << (2*order_);
  
  long ncap_ = (npface_-nside_)<<1;
  
  long nl4 = 4*nside_;
  long jr = (jrll[face_num]*nside_) - ix - iy  - 1;

  long nr, kshift, n_before;
  if (jr<nside_)
    {
      nr = jr;
      n_before = 2*nr*(nr-1);
      kshift = 0;
    }
  else if (jr > 3*nside_)
    {
      nr = nl4-jr;
      n_before = npix_ - 2*(nr+1)*nr;
      kshift = 0;
    }
  else
    {
      nr = nside_;
      n_before = ncap_ + (jr-nside_)*nl4;
      kshift = (jr-nside_)&1;
    }

  long jp = (jpll[face_num]*nr + ix - iy + 1 + kshift) / 2;
  if (jp>nl4)
    jp-=nl4;
  else
    if (jp<1) jp+=nl4;

  long pix = n_before + jp - 1;
  assert(pix >= 0 && pix < npix_);
  
  return pix;
}

void nest2xyf(long pix, long *ix, long *iy, long *face_num, long order_)
{
  if(HEALPIX_TOOLS_INIT)
    {
      tablefiller();
      HEALPIX_TOOLS_INIT = 0;
    }
  
  long npix_ = 1;
  npix_ = 12*(npix_ << (2*order_));
  assert(pix >= 0 && pix < npix_);
  
  long npface_ = 1;
  npface_ = npface_ << (2*order_);
  
  *face_num = pix>>(2*order_);
  pix &= (npface_-1);
  long raw = ((pix&0x555500000000ull)>>16) 
    | ((pix&0x5555000000000000ull)>>31)
    | (pix&0x5555)
    | ((pix&0x55550000)>>15);
  *ix =  ctab[raw&0xff]
    | (ctab[(raw>>8)&0xff]<<4)
    | (ctab[(raw>>16)&0xff]<<16)
    | (ctab[(raw>>24)&0xff]<<20);
  pix >>= 1;
  raw = ((pix&0x555500000000ull)>>16) 
    | ((pix&0x5555000000000000ull)>>31)
    | (pix&0x5555)
    | ((pix&0x55550000)>>15);
  *iy =  ctab[raw&0xff]
    | (ctab[(raw>>8)&0xff]<<4)
    | (ctab[(raw>>16)&0xff]<<16)
    | (ctab[(raw>>24)&0xff]<<20);
}


void ring2xyf(long pix, long *ix, long *iy, long *face_num, long order_)
{
  long iring, iphi, kshift, nr;
  
  long npix_ = 1;
  npix_ = 12*(npix_ << (2*order_));
  assert(pix >= 0 && pix < npix_);
  
  long nside_ = 1;
  nside_ = nside_ << order_;
    
  long npface_ = 1;
  npface_ = npface_ << (2*order_);
  
  long ncap_ = (npface_-nside_)<<1;
  
  long nl2 = 2*nside_;
  
  if (pix<ncap_) // North Polar cap
    {
      iring = (long) (0.5*(1+isqrt(1+2*pix))); //counted from North pole
      iphi  = (pix+1) - 2*iring*(iring-1);
      kshift = 0;
      nr = iring;
      *face_num=0;
      long tmp = iphi-1;
      if (tmp>=(2*iring))
	{
	  *face_num=2;
	  tmp-=2*iring;
	}
      if (tmp>=iring) (*face_num) = (*face_num) + 1;
    }
  else if (pix<(npix_-ncap_)) // Equatorial region
    {
      long ip = pix - ncap_;
      if (order_>=0)
	{
	  iring = (ip>>(order_+2)) + nside_; // counted from North pole
	  iphi  = (ip&(4*nside_-1)) + 1;
	}
      else
	{
	  iring = (ip/(4*nside_)) + nside_; // counted from North pole
	  iphi  = (ip%(4*nside_)) + 1;
	}
      kshift = (iring+nside_)&1;
      nr = nside_;
      long ire = iring-nside_+1;
      long irm = nl2+2-ire;
      long ifm, ifp;
      if (order_>=0)
	{
	  ifm = (iphi - ire/2 + nside_ -1) >> order_;
	  ifp = (iphi - irm/2 + nside_ -1) >> order_;
	}
      else
	{
	  ifm = (iphi - ire/2 + nside_ -1) / nside_;
	  ifp = (iphi - irm/2 + nside_ -1) / nside_;
	}
      if (ifp == ifm) // faces 4 to 7
	*face_num = (ifp==4) ? 4 : ifp+4;
      else if (ifp<ifm) // (half-)faces 0 to 3
	*face_num = ifp;
      else // (half-)faces 8 to 11
	*face_num = ifm + 8;
    }
  else // South Polar cap
    {
      long ip = npix_ - pix;
      iring = (long) (0.5*(1+isqrt(2*ip-1))); //counted from South pole
      iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1));
      kshift = 0;
      nr = iring;
      iring = 2*nl2-iring;
      *face_num=8;
      long tmp = iphi-1;
      if (tmp>=(2*nr))
	{
	  *face_num=10;
	  tmp-=2*nr;
	}
      if (tmp>=nr) (*face_num) = (*face_num) + 1;
    }

  long irt = iring - (jrll[*face_num]*nside_) + 1;
  long ipt = 2*iphi- jpll[*face_num]*nr - kshift -1;
  if (ipt>=nl2) ipt-=8*nside_;

  *ix =  (ipt-irt) >>1;
  *iy =(-(ipt+irt))>>1;

}

void ang2pix_ring(long nside, double theta, double phi, long *pix)
{
  int order_, temp;

  temp = nside;
  order_ = 0;
  while(temp >>= 1) order_++;

  long npix_ = 1;
  npix_ = 12*(npix_ << (2*order_));
  
  long nside_ = 1;
  nside_ = nside_ << order_;
  
  long npface_ = 1;
  npface_ = npface_ << (2*order_);
  
  long ncap_ = (npface_-nside_)<<1;
  
  double z = cos(theta);
  double za = fabs(z);
  double tt = phi;
  long tt_long = (floor(tt/2/M_PI));
  tt = tt - ((double) (tt_long))*2*M_PI;
  tt *= M_2_PI; // in [0,4)

  if(za <= 2.0/3.0) // Equatorial region
    {
      double temp1 = nside_*(0.5+tt);
      double temp2 = nside_*z*0.75;
      long jp = (long) (temp1-temp2); // index of  ascending edge line
      long jm = (long) (temp1+temp2); // index of descending edge line
      
      // ring number counted from z=2/3
      long ir = nside_ + 1 + jp - jm; // in {1,2n+1}
      long kshift = 1-(ir&1); // kshift=1 if ir even, 0 otherwise
      
      long ip = (jp+jm-nside_+kshift+1)/2; // in {0,4n-1}
      ip = imodulo(ip,4*nside_);
      
      *pix = ncap_ + (ir-1)*4*nside_ + ip;
    }
  else  // North & South polar caps
    {
      double tp = tt- ((long) (tt));
      double tmp = nside_*sqrt(3*(1-za));
      
      long jp = (long) (tp*tmp); // increasing edge line index
      long jm = (long) ((1.0-tp)*tmp); // decreasing edge line index
      
      long ir = jp+jm+1; // ring number counted from the closest pole
      long ip = (long) (tt*ir); // in {0,4*ir-1}
      ip = imodulo(ip,4*ir);
      
      if(z>0)
	*pix = 2*ir*(ir-1) + ip;
      else
	*pix = npix_ - 2*ir*(ir+1) + ip;
    }
}
