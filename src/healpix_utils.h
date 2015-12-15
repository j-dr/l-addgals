#ifndef HEALPIXUTILS_H
#define HEALPIXUTILS_H
#include <vector>

long nest2peano(long pix, long order_);
void ring2peanoindex(long pix, long order1_, long order2_, std::vector<long> &pidx);
long ring2nest(long pix, long order_);
long nest2ring(long pix, long order_);
void higher_nest(long pix, long order1_, long order2_, long* hopix);
void nest2xyf(long pix, long *ix, long *iy, long *face_num, long order_);
long xyf2nest(long ix, long iy, long face_num, long order_);
void ring2xyf(long pix, long *ix, long *iy, long *face_num, long order_);
long xyf2ring(long ix, long iy, long face_num, long order_);
void tablefiller(void);


#endif
