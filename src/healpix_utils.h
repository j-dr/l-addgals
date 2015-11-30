#ifndef HEALPIXUTILS_H
#define HEALPIXUTILS_H

long nest2peano(long pix, long order_);
long ring2nest(long pix, long order_);
long nest2ring(long pix, long order_);
long higher_nest(long pix, long order1_, long order2_);
void nest2xyf(long pix, long *ix, long *iy, long *face_num, long order_);
long xyf2nest(long ix, long iy, long face_num, long order_);
void ring2xyf(long pix, long *ix, long *iy, long *face_num, long order_);
long xyf2ring(long ix, long iy, long face_num, long order_)


#endif
