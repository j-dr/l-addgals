#include <iostream>
#include <fstream>
#include <vector>
#include <ANN/ANN.h>			// ANN declarations
#include "nr.h"
#include "galaxy.h"
#include "singleton.h"
#include "choose.h"
#include <math.h>

int findCloseGalaxies2(vector <GalSED> &v, float mag, float den, float ThisZ, int ThisBCGs);

/*

#include <ANN/ANNx.h>			// all ANN includes
#include "/home/risa/code/ann_0.2/src/kd_tree.h"			// kd-tree declarations
#include "/home/risa/code/ann_0.2/src/kd_search.h"			// kd-search declarations
#include "/home/risa/code/ann_0.2/src/kd_split.h"			// kd-tree splitting rules
#include "/home/risa/code/ann_0.2/src/kd_util.h"			// kd-tree utilities
//#include <ANN/ANNperf.h>		// performance evaluation

static int  		IDX_TRIVIAL[] = {0};	// trivial point index
ANNkd_leaf		*KD_TRIVIAL = NULL;	// trivial leaf node


ANNpoint annAllocPt(int dim, ANNcoord c)	// allocate 1 point
{
    ANNpoint p = new ANNcoord[dim];
    for (int i = 0; i < dim; i++) p[i] = c;
    return p;
}
   
ANNpointArray annAllocPts(int n, int dim)	// allocate n pts in dim
{
    ANNpointArray pa = new ANNpoint[n];		// allocate points
    ANNpoint	  p  = new ANNcoord[n*dim];	// allocate space for coords
    for (int i = 0; i < n; i++) {
	pa[i] = &(p[i*dim]);
    }
    return pa;
}

void annDeallocPt(ANNpoint &p)			// deallocate 1 point
{
    delete [] p;
    p = NULL;
}
   
void annDeallocPts(ANNpointArray &pa)		// deallocate points
{
    delete [] pa[0];				// dealloc coordinate storage
    delete [] pa;				// dealloc points
    pa = NULL;
}
  

int             ANNmaxPtsVisited = 0;   // maximum number of pts visited
int             ANNptsVisited;          // number of pts visited in search



ANNpoint annCopyPt(int dim, ANNpoint source)    // copy point
{
    ANNpoint p = new ANNcoord[dim];
    for (int i = 0; i < dim; i++) p[i] = source[i];
    return p;
}

                                                // assign one rect to another
void annAssignRect(int dim, ANNorthRect &dest, const ANNorthRect &source)
{
    for (int i = 0; i < dim; i++) {
        dest.lo[i] = source.lo[i];
        dest.hi[i] = source.hi[i];
    }
}



                                                // is point inside rectangle?
ANNbool ANNorthRect::inside(int dim, ANNpoint p)
{
    for (int i = 0; i < dim; i++) {
        if (p[i] < lo[i] || p[i] > hi[i]) return ANNfalse;
    }
    return ANNtrue;
}

//----------------------------------------------------------------------
//  Error handler
//----------------------------------------------------------------------

void annError(char *msg, ANNerr level)
{
    if (level == ANNabort) {
        cerr << "ANN: ERROR------->" << msg << "<-------------ERROR\n";
        exit(1);
    }
    else {
        cerr << "ANN: WARNING----->" << msg << "<-------------WARNING\n";
    }
}



//----------------------------------------------------------------------
//  Approximate nearest neighbor searching by kd-tree search
//	The kd-tree is searched for an approximate nearest neighbor.
//	The point is returned through one of the arguments, and the
//	distance returned is the squared distance to this point.
//
//	The method used for searching the kd-tree is an approximate
//	adaptation of the search algorithm described by Friedman,
//	Bentley, and Finkel, ``An algorithm for finding best matches
//	in logarithmic expected time,'' ACM Transactions on Mathematical
//	Software, 3(3):209-226, 1977).
//
//	The algorithm operates recursively.  When first encountering a
//	node of the kd-tree we first visit the child which is closest to
//	the query point.  On return, we decide whether we want to visit
//	the other child.  If the box containing the other child exceeds
//	1/(1+eps) times the current best distance, then we skip it (since
//	any point found in this child cannot be closer to the query point
//	by more than this factor.)  Otherwise, we visit it recursively.
//	The distance between a box and the query point is computed exactly
//	(not approximated as is often done in kd-tree), using incremental
//	distance updates, as described by Arya and Mount in ``Algorithms
//	for fast vector quantization,'' Proc.  of DCC '93: Data Compression
//	Conference, eds. J. A. Storer and M. Cohn, IEEE Press, 1993,
//	381-390.
//
//	The main entry points is annkSearch() which sets things up and
//	then call the recursive routine ann_search().  This is a recursive
//	routine which performs the processing for one node in the kd-tree.
//	There are two versions of this virtual procedure, one for splitting
//	nodes and one for leaves.  When a splitting node is visited, we
//	determine which child to visit first (the closer one), and visit
//	the other child on return.  When a leaf is visited, we compute
//	the distances to the points in the buckets, and update information
//	on the closest points.
//
//	Some trickery is used to incrementally update the distance from
//	a kd-tree rectangle to the query point.  This comes about from
//	the fact that which each successive split, only one component
//	(along the dimension that is split) of the squared distance to
//	the child rectangle is different from the squared distance to
//	the parent rectangle.
//----------------------------------------------------------------------

//----------------------------------------------------------------------
//	To keep argument lists short, a number of global variables
//	are maintained which are common to all the recursive calls.
//	These are given below.
//----------------------------------------------------------------------

int		ANNkdDim;		// dimension of space
ANNpoint	ANNkdQ;			// query point
double		ANNkdMaxErr;		// max tolerable squared error
ANNpointArray	ANNkdPts;		// the points
ANNmin_k	*ANNkdPointMK;		// set of k closest points

//----------------------------------------------------------------------
//  annkSearch - search for the k nearest neighbors
//----------------------------------------------------------------------

void ANNkd_tree::SkeletonTree(          // construct skeleton tree
        int n,                          // number of points
        int dd,                         // dimension
        int bs)                         // bucket size
{
    dim = dd;                           // initialize basic elements
    n_pts = n;
    bkt_size = bs;
    root = NULL;                        // no associated tree yet
    pts = NULL;                         // no associated points yet

    pidx = new int[n];                  // allocate space for point indices
    for (int i = 0; i < n; i++) {
        pidx[i] = i;                    // initially identity
    }
    bnd_box_lo = bnd_box_hi = NULL;     // bounding box is nonexistent
    if (KD_TRIVIAL == NULL)             // no trivial leaf node yet?
        KD_TRIVIAL = new ANNkd_leaf(0, IDX_TRIVIAL);    // allocate it
}



ANNkd_tree::ANNkd_tree(                 // basic constructor
        int n,                          // number of points
        int dd,                         // dimension
        int bs)                         // bucket size
{  SkeletonTree(n, dd, bs);            // construct skeleton tree
}

//----------------------------------------------------------------------
// kd-tree constructor
//      This is the main constructor for kd-trees given a set of points.
//      It first builds a skeleton tree, then computes the bounding box
//      of the data points, and then invokes rkd_tree() to actually
//      build the tree, passing it the appropriate splitting routine.
//----------------------------------------------------------------------

ANNkd_tree::ANNkd_tree(                 // construct from point array
    ANNpointArray       pa,             // point array (with at least n pts)
    int                 n,              // number of points
    int                 dd,             // dimension
    int                 bs,             // bucket size
    ANNsplitRule        split)          // splitting method
{
    SkeletonTree(n, dd, bs);            // set up the basic stuff
    pts = pa;                           // where the points are
    if (n == 0) return;                 // no points--no sweat

    ANNorthRect bnd_box(dd);            // bounding box for points
    annEnclRect(pa, pidx, n, dd, bnd_box);// construct bounding rectangle
                                        // copy to tree structure
    bnd_box_lo = annCopyPt(dd, bnd_box.lo);
    bnd_box_hi = annCopyPt(dd, bnd_box.hi);

    switch (split) {                    // build by rule
    case ANN_KD_STD:                    // standard kd-splitting rule
        root = rkd_tree(pa, pidx, n, dd, bs, bnd_box, kd_split);
        break;
    case ANN_KD_MIDPT:                  // midpoint split
        root = rkd_tree(pa, pidx, n, dd, bs, bnd_box, midpt_split);
        break;
    case ANN_KD_FAIR:                   // fair split
        root = rkd_tree(pa, pidx, n, dd, bs, bnd_box, fair_split);
        break;
    case ANN_KD_SUGGEST:                // best (in our opinion)
    case ANN_KD_SL_MIDPT:               // sliding midpoint split
        root = rkd_tree(pa, pidx, n, dd, bs, bnd_box, sl_midpt_split);
        break;
    case ANN_KD_SL_FAIR:                // sliding fair split
        root = rkd_tree(pa, pidx, n, dd, bs, bnd_box, sl_fair_split);
        break;
    default:
        annError("Illegal splitting method", ANNabort);
    }
}






void ANNkd_tree::annkSearch(
    ANNpoint		q,		// the query point
    int			k,		// number of near neighbors to return
    ANNidxArray		nn_idx,		// nearest neighbor indices (returned)
    ANNdistArray	dd,		// the approximate nearest neighbor
    double		eps)		// the error bound
{

    ANNkdDim = dim;			// copy arguments to static equivs
    ANNkdQ = q;
    ANNkdPts = pts;
    ANNptsVisited = 0;			// initialize count of points visited

    if (k > n_pts) {			// too many near neighbors?
      cerr<<k<<" "<<n_pts<<endl;
      annError("Requesting more near neighbors than data points", ANNabort);
    }

    ANNkdMaxErr = ANN_POW(1.0 + eps);
    FLOP(2)				// increment floating op count

    ANNkdPointMK = new ANNmin_k(k);	// create set for closest k points
					// search starting at the root
    root->ann_search(annBoxDistance(q, bnd_box_lo, bnd_box_hi, dim));

    for (int i = 0; i < k; i++) {	// extract the k-th closest points
	dd[i] = ANNkdPointMK->ith_smallest_key(i);
	nn_idx[i] = ANNkdPointMK->ith_smallest_info(i);
    }
    delete ANNkdPointMK;		// deallocate closest point set
}

//----------------------------------------------------------------------
//  kd_split::ann_search - search a splitting node
//----------------------------------------------------------------------

void ANNkd_split::ann_search(ANNdist box_dist)
{
					// check dist calc termination condition
    if (ANNmaxPtsVisited && ANNptsVisited > ANNmaxPtsVisited) return;

					// distance to cutting plane
    ANNcoord cut_diff = ANNkdQ[cut_dim] - cut_val;

    if (cut_diff < 0) {			// left of cutting plane
	child[LO]->ann_search(box_dist);// visit closer child first

	ANNcoord box_diff = cd_bnds[LO] - ANNkdQ[cut_dim];
	if (box_diff < 0)		// within bounds - ignore
	    box_diff = 0;
					// distance to further box
	box_dist = (ANNdist) ANN_SUM(box_dist,
		ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

					// visit further child if close enough
        if (box_dist * ANNkdMaxErr < ANNkdPointMK->max_key())
            child[HI]->ann_search(box_dist);

    }
    else {				// right of cutting plane
	child[HI]->ann_search(box_dist);// visit closer child first

	ANNcoord box_diff = ANNkdQ[cut_dim] - cd_bnds[HI];
	if (box_diff < 0)		// within bounds - ignore
	    box_diff = 0;
					// distance to further box
	box_dist = (ANNdist) ANN_SUM(box_dist,
		ANN_DIFF(ANN_POW(box_diff), ANN_POW(cut_diff)));

					// visit further child if close enough
        if (box_dist * ANNkdMaxErr < ANNkdPointMK->max_key())
            child[LO]->ann_search(box_dist);

    }
    FLOP(10)				// increment floating ops
    SPL(1)				// one more splitting node visited
}

//----------------------------------------------------------------------
//  kd_leaf::ann_search - search points in a leaf node
//	Note: The unreadability of this code is the result of
//	some fine tuning to replace indexing by pointer operations.
//----------------------------------------------------------------------

void ANNkd_leaf::ann_search(ANNdist box_dist)
{
    register ANNdist dist;		// distance to data point
    register ANNcoord* pp;		// data coordinate pointer
    register ANNcoord* qq;		// query coordinate pointer
    register ANNdist min_dist;		// distance to k-th closest point
    register ANNcoord t;
    register int d;

    min_dist = ANNkdPointMK->max_key();	// k-th smallest distance so far

    for (int i = 0; i < n_pts; i++) {	// check points in bucket

	pp = ANNkdPts[bkt[i]];		// first coord of next data point
	qq = ANNkdQ; 	    		// first coord of query point
	dist = 0;

	for(d = 0; d < ANNkdDim; d++) {
	    COORD(1)			// one more coordinate hit
	    FLOP(4)			// increment floating ops

	    t = *(qq++) - *(pp++);	// compute length and adv coordinate
					// exceeds dist to k-th smallest?
	    if( (dist = ANN_SUM(dist, ANN_POW(t))) > min_dist) {
		break;
	    }
	}

	if (d >= ANNkdDim &&			// among the k best?
	   (ANN_ALLOW_SELF_MATCH || dist!=0)) {	// and no self-match problem
						// add it to the list
	    ANNkdPointMK->insert(dist, bkt[i]);
	    min_dist = ANNkdPointMK->max_key();
	}
    }
    LEAF(1)				// one more leaf node visited
    PTS(n_pts)				// increment points visited
    ANNptsVisited += n_pts;		// increment number of points visited
}

//----------------------------------------------------------------------
//	File:		kd_util.cc
//	Programmer:	Sunil Arya and David Mount
//	Last modified:	03/04/98 (Release 0.1)
//	Description:	Common utilities for kd-trees
//
//	The following routines are utility functions for manipulating
//	points sets, used in determining splitting planes for kd-tree
//	construction.
//----------------------------------------------------------------------
// Copyright (c) 1997-1998 University of Maryland and Sunil Arya and David
// Mount.  All Rights Reserved.
// 
// This software and related documentation is part of the 
// Approximate Nearest Neighbor Library (ANN).
// 
// Permission to use, copy, and distribute this software and its 
// documentation is hereby granted free of charge, provided that 
// (1) it is not a component of a commercial product, and 
// (2) this notice appears in all copies of the software and
//     related documentation. 
// 
// The University of Maryland (U.M.) and the authors make no representations
// about the suitability or fitness of this software for any purpose.  It is
// provided "as is" without express or implied warranty.
//----------------------------------------------------------------------


//----------------------------------------------------------------------
//	NOTE: Virtually all point indexing is done through an
//	index (i.e. permutation) array pidx.  Consequently,
//	a reference to the d-th coordinate of the i-th point
//	is pa[pidx[i]][d].  The macro PA(i,d) is a shorthand
//	for this.
//----------------------------------------------------------------------
					// standard 2-d indirect indexing
#define PA(i,d) 	(pa[pidx[(i)]][(d)])
					// accessing a single point
#define PP(i)	 	(pa[pidx[(i)]])

//----------------------------------------------------------------------
//  annAspectRatio
//	Compute the aspect ratio (ratio of longest to shortest side)
//	of a rectangle.
//----------------------------------------------------------------------

double annAspectRatio(
    int			dim,		// dimension
    const ANNorthRect	&bnd_box)	// bounding cube
{
    ANNcoord length = bnd_box.hi[0] - bnd_box.lo[0];
    ANNcoord min_length = length;		// min side length
    ANNcoord max_length = length;		// max side length
    for (int d = 0; d < dim; d++) {
	length = bnd_box.hi[d] - bnd_box.lo[d];
	if (length < min_length) min_length = length;
	if (length > max_length) max_length = length;
    }
    return max_length/min_length;
}

//----------------------------------------------------------------------
//  annEnclRect, annEnclCube
//	These utilities compute the smallest rectangle and cube enclosing
//	a set of points, respectively.
//----------------------------------------------------------------------

void annEnclRect(
    ANNpointArray	pa,		// point array
    ANNidxArray		pidx,		// point indices
    int			n,		// number of points
    int			dim,		// dimension
    ANNorthRect		&bnds)		// bounding cube (returned)
{
    for (int d = 0; d < dim; d++) {	// find smallest enclosing rectangle
	ANNcoord lo_bnd = PA(0,d);	// lower bound on dimension d
	ANNcoord hi_bnd = PA(0,d);	// upper bound on dimension d
        for (int i = 0; i < n; i++) {
	    if (PA(i,d) < lo_bnd) lo_bnd = PA(i,d);
	    else if (PA(i,d) > hi_bnd) hi_bnd = PA(i,d);
	}
	bnds.lo[d] = lo_bnd;
	bnds.hi[d] = hi_bnd;
    }
}

void annEnclCube(			// compute smallest enclosing cube
    ANNpointArray	pa,		// point array
    ANNidxArray		pidx,		// point indices
    int			n,		// number of points
    int			dim,		// dimension
    ANNorthRect		&bnds)		// bounding cube (returned)
{
    int d;
					// compute smallest enclosing rect
    annEnclRect(pa, pidx, n, dim, bnds);

    ANNcoord max_len = 0;		// max length of any side
    for (d = 0; d < dim; d++) {		// determine max side length
	ANNcoord len = bnds.hi[d] - bnds.lo[d];
	if (len > max_len) {		// update max_len if longest
	    max_len = len;
	}
    }
    for (d = 0; d < dim; d++) {		// grow sides to match max
	ANNcoord len = bnds.hi[d] - bnds.lo[d];
	ANNcoord half_diff = (max_len - len) / 2;
	bnds.lo[d] -= half_diff;
	bnds.hi[d] += half_diff;
    }
}

//----------------------------------------------------------------------
//  annBoxDistance - utility routine which computes distance from point to
//	box (Note: most distances to boxes are computed using incremental
//	distance updates, not this function.)
//----------------------------------------------------------------------

ANNdist annBoxDistance(		// compute distance from point to box
    const ANNpoint	q,		// the point
    const ANNpoint	lo,		// low point of box
    const ANNpoint	hi,		// high point of box
    int			dim)		// dimension of space
{
    register ANNdist dist = 0.0;	// sum of squared distances
    register ANNdist t;

    for (register int d = 0; d < dim; d++) {
	if (q[d] < lo[d]) {		// q is left of box
	    t = ANNdist(lo[d]) - ANNdist(q[d]);
	    dist = ANN_SUM(dist, ANN_POW(t));
	}
	else if (q[d] > hi[d]) {	// q is right of box
	    t = ANNdist(q[d]) - ANNdist(hi[d]);
	    dist = ANN_SUM(dist, ANN_POW(t));
	}
    }
    FLOP(4*dim)				// increment floating op count

    return dist;
}

//----------------------------------------------------------------------
//  annSpread - find spread along given dimension
//  annMinMax - find min and max coordinates along given dimension
//  annMaxSpread - find dimension of max spread
//----------------------------------------------------------------------

ANNcoord annSpread(		// compute point spread along dimension
    ANNpointArray	pa,		// point array
    ANNidxArray		pidx,		// point indices
    int			n,		// number of points
    int			d)		// dimension to check
{
    ANNcoord min = PA(0,d);		// compute max and min coords
    ANNcoord max = PA(0,d);
    for (int i = 1; i < n; i++) {
	ANNcoord c = PA(i,d);
	if (c < min) min = c;
	else if (c > max) max = c;
    }
    return (max - min);			// total spread is difference
}

void annMinMax(			// compute min and max coordinates along dim
    ANNpointArray	pa,		// point array
    ANNidxArray		pidx,		// point indices
    int			n,		// number of points
    int			d,		// dimension to check
    ANNcoord		&min,		// minimum value (returned)
    ANNcoord		&max)		// maximum value (returned)
{
    min = PA(0,d);			// compute max and min coords
    max = PA(0,d);
    for (int i = 1; i < n; i++) {
	ANNcoord c = PA(i,d);
	if (c < min) min = c;
	else if (c > max) max = c;
    }
}

int annMaxSpread(			// compute dimension of max spread
    ANNpointArray	pa,		// point array
    ANNidxArray		pidx,		// point indices
    int			n,		// number of points
    int			dim)		// dimension of space
{
    int max_dim = 0;			// dimension of max spread
    ANNcoord max_spr = 0;		// amount of max spread

    if (n == 0) return max_dim;		// no points, who cares?

    for (int d = 0; d < dim; d++) {	// compute spread along each dim
	ANNcoord spr = annSpread(pa, pidx, n, d);
	if (spr > max_spr) {		// bigger than current max
	    max_spr = spr;
	    max_dim = d;
	}
    }
    return max_dim;
}

//----------------------------------------------------------------------
//  annMedianSplit - split point array about its median
//	Splits a subarray of points pa[0..n] about an element of given
//	rank (median: n_lo = n/2) with respect to dimension d.  It places
//	the element of rank n_lo-1 correctly (because our splitting rule
//	takes the mean of these two).  On exit, the array is permuted so
//	that:
//
//	pa[0..n_lo-2][d] <= pa[n_lo-1][d] <= pa[n_lo][d] <= pa[n_lo+1..n-1][d].
//
//	The mean of pa[n_lo-1][d] and pa[n_lo][d] is returned as the
//	splitting value.
//
//	All indexing is done indirectly through the index array pidx.
//
//	This function uses the well known selection algorithm due to
//	C.A.R. Hoare.
//----------------------------------------------------------------------
					// swap two points in pa array
#define PASWAP(a,b) { int tmp = pidx[a];\
		    pidx[a] = pidx[b];\
		    pidx[b] = tmp; }

void annMedianSplit(
    ANNpointArray	pa,		// points to split
    ANNidxArray		pidx,		// point indices
    int			n,		// number of points
    int			d,		// dimension along which to split
    ANNcoord		&cv,		// cutting value
    int			n_lo)		// split into n_lo and n-n_lo
{
    int l = 0;				// left end of current subarray
    int r = n-1;			// right end of current subarray
    while (l < r) {
	register int i = (r+l)/2;	// select middle as pivot
	register int k;

	if (PA(i,d) > PA(r,d))		// make sure last > pivot
	    PASWAP(i,r)
	PASWAP(l,i);			// move pivot to first position

	ANNcoord c = PA(l,d);		// pivot value
	i = l;
	k = r;
	for(;;) {			// pivot about c
	    while (PA(++i,d) < c) ;
	    while (PA(--k,d) > c) ;
	    if (i < k) PASWAP(i,k) else break;
	}
	PASWAP(l,k);			// pivot winds up in location k

	if (k > n_lo)      r = k-1;	// recurse on proper subarray
	else if (k < n_lo) l = k+1;
	else break;			// got the median exactly
    }
    if (n_lo > 0) {			// search for next smaller item
	ANNcoord c = PA(0,d);		// candidate for max
	int k = 0;			// candidate's index
	for (int i = 1; i < n_lo; i++) {
	    if (PA(i,d) > c) {
		c = PA(i,d);
		k = i;
	    }
	}
	PASWAP(n_lo-1, k);		// max among pa[0..n_lo-1] to pa[n_lo-1]
    }
					// cut value is midpoint value
    cv = (PA(n_lo-1,d) + PA(n_lo,d))/2.0;
}

//----------------------------------------------------------------------
//  annPlaneSplit - split point array about a cutting plane
//	Split the points in an array about a given plane along a
//	given cutting dimension.  On exit, br1 and br2 are set so
//	that:
//	
//		pa[ 0 ..br1-1] <  cv
//		pa[br1..br2-1] == cv
//		pa[br2.. n -1] >  cv
//
//	All indexing is done indirectly through the index array pidx.
//
//----------------------------------------------------------------------

void annPlaneSplit(		// split points by a plane
    ANNpointArray	pa,		// points to split
    ANNidxArray		pidx,		// point indices
    int			n,		// number of points
    int			d,		// dimension along which to split
    ANNcoord		cv,		// cutting value
    int			&br1,		// first break (values < cv)
    int			&br2)		// second break (values == cv)
{
    int l = 0;
    int r = n-1;
    for(;;) {				// partition pa[0..n-1] about cv
	while (l < n && PA(l,d) < cv) l++;
	while (r >= 0 && PA(r,d) >= cv) r--;
	if (l > r) break;
	PASWAP(l,r);
	l++; r--;
    }
    br1 = l;			// now: pa[0..br1-1] < cv <= pa[br1..n-1]
    r = n-1;
    for(;;) {				// partition pa[br1..n-1] about cv
	while (l < n && PA(l,d) <= cv) l++;
	while (r >= br1 && PA(r,d) > cv) r--;
	if (l > r) break;
	PASWAP(l,r);
	l++; r--;
    }
    br2 = l;			// now: pa[br1..br2-1] == cv < pa[br2..n-1]
}


//----------------------------------------------------------------------
//  annBoxSplit - split point array about a orthogonal rectangle
//	Split the points in an array about a given orthogonal
//	rectangle.  On exit, n_in is set to the number of points
//	that are inside (or on the boundary of) the rectangle.
//
//	All indexing is done indirectly through the index array pidx.
//
//----------------------------------------------------------------------

void annBoxSplit(		// split points by a box
    ANNpointArray	pa,		// points to split
    ANNidxArray		pidx,		// point indices
    int			n,		// number of points
    int			dim,		// dimension of space
    ANNorthRect		&box,		// the box
    int			&n_in)		// number of points inside (returned)
{
    int l = 0;
    int r = n-1;
    for(;;) {				// partition pa[0..n-1] about box
	while (l < n && box.inside(dim, PP(l))) l++;
	while (r >= 0 && !box.inside(dim, PP(r))) r--;
	if (l > r) break;
	PASWAP(l,r);
	l++; r--;
    }
    n_in = l;			// now: pa[0..n_in-1] inside and rest outside
}

//----------------------------------------------------------------------
//  annSplitBalance - compute balance factor for a given plane split
//	Balance factor is defined as the number of points lying
//	below the splitting value minus n/2 (median).  Thus, a
//	median split has balance 0, left of this is negative and
//	right of this is positive.  (The points are unchanged.)
//----------------------------------------------------------------------

int annSplitBalance(		// determine balance factor of a split
    ANNpointArray	pa,		// points to split
    ANNidxArray		pidx,		// point indices
    int			n,		// number of points
    int			d,		// dimension along which to split
    ANNcoord		cv)		// cutting value
{
    int n_lo = 0;
    for(int i = 0; i < n; i++) {	// count number less than cv
	if (PA(i,d) < cv) n_lo++;
    }
    return n_lo - n/2;
}

//----------------------------------------------------------------------
//  annBox2Bnds - convert bounding box to list of bounds
//	Given two boxes, an inner box enclosed within a bounding
//	box, this routine determines all the sides for which the
//	inner box is strictly contained with the bounding box,
//	and adds an appropriate entry to a list of bounds.  Then
//	we allocate storage for the final list of bounds, and return
//	the resulting list and its size.
//----------------------------------------------------------------------

void annBox2Bnds(			// convert inner box to bounds
    const ANNorthRect	&inner_box,	// inner box
    const ANNorthRect	&bnd_box,	// enclosing box
    int			dim,		// dimension of space
    int			&n_bnds,	// number of bounds (returned)
    ANNorthHSArray	&bnds)		// bounds array (returned)
{
    int i;
    n_bnds = 0;					// count number of bounds
    for (i = 0; i < dim; i++) {
	if (inner_box.lo[i] > bnd_box.lo[i])	// low bound is inside
		n_bnds++;
	if (inner_box.hi[i] < bnd_box.hi[i])	// high bound is inside
		n_bnds++;
    }

    bnds = new ANNorthHalfSpace[n_bnds];	// allocate appropriate size

    int j = 0;
    for (i = 0; i < dim; i++) {			// fill the array
	if (inner_box.lo[i] > bnd_box.lo[i]) {
		bnds[j].cd = i;
		bnds[j].cv = inner_box.lo[i];
		bnds[j].sd = +1;
		j++;
	}
	if (inner_box.hi[i] < bnd_box.hi[i]) {
		bnds[j].cd = i;
		bnds[j].cv = inner_box.hi[i];
		bnds[j].sd = -1;
		j++;
	}
    }
}

//----------------------------------------------------------------------
//  annBnds2Box - convert list of bounds to bounding box
//	Given an enclosing box and a list of bounds, this routine
//	computes the corresponding inner box.  It is assumed that
//	the box points have been allocated already.
//----------------------------------------------------------------------

void annBnds2Box(
    const ANNorthRect	&bnd_box,	// enclosing box
    int			dim,		// dimension of space
    int			n_bnds,		// number of bounds
    ANNorthHSArray	bnds,		// bounds array
    ANNorthRect		&inner_box)	// inner box (returned)
{
    annAssignRect(dim, inner_box, bnd_box);	// copy bounding box to inner

    for (int i = 0; i < n_bnds; i++) {
	bnds[i].project(inner_box.lo);		// project each endpoint
	bnds[i].project(inner_box.hi);
    }
}

*/


vector <GalSED> ReadSED();
/*
// new 2d version                                                             
// this calculates the proper ra-dec distance                                 
vector <float> GetNeighborDistAM(vector <Galaxy *> galaxies){
  vector <float> nndist(galaxies.size());
  int                 n_pts = galaxies.size();                  // actual number of data points                                                            
  ANNpointArray       data_pts;               // data points                  
 ANNpoint            query_pt;               // query point                  
 ANNidxArray         nn_idx;                 // near neighbor indices        
 ANNdistArray        dists;                  // near neighbor distances      
 ANNkd_tree          *the_tree;              // search structure             
 double eps=0;
 int dim = 2;
 int k_want = 11;  //counts itself, so this is really the 10th nearest neighbor                                                                           
  if(dens_measure == FIFTH)
    k_want = 5;
 int k = 100;
 query_pt = annAllocPt(dim);                 // allocate query point         
 data_pts = annAllocPts(n_pts, dim);         // allocate data points         
 nn_idx = new ANNidx[k];                     // allocate near neigh indices  
 dists = new ANNdist[k];                     // allocate near neighbor dists

 int jmax = 0;

 //ofstream outpts("pts.out");                                               
 for(int i=0;i<n_pts;i++){
   //Particle *p = galaxies[i]->P();                                         
   data_pts[i][0] = galaxies[i]->Ra();
   data_pts[i][1] = galaxies[i]->Dec();
   //data_pts[i][0] = p->X();                                                
   //data_pts[i][1] = p->Y();                                                
   //outpts<<data_pts[i][0]<<" "<<data_pts[i][1]<<endl;                        
 }

 //  ofstream outall("all.out");                                             
 the_tree = new ANNkd_tree(                    // build search structure     
			   data_pts,                   // the data points    
			   n_pts,                      // number of points   
			   dim);                       // dimension of space

 ofstream out("dd2.out");

 for(int i=0;i<n_pts;i++){
   if(i%20000==0) {cout<<i<<endl; system("date");}
   bool distfound = false;

   ANNpoint query_pt = data_pts[i];
   float dist=0;
   k=100;
   while(!distfound){
     int nn = 0;
     //if(k>=800) cout<<"searching..."<<k<<" ";                              
     the_tree->annkSearch(                     // search                     
			  query_pt,            // query point                
			  k,                   // number of near neighbors   
			  nn_idx,              // nearest neighbors (returned)                                                                            
       dists,               // distance (returned)        
       eps);                // error bound                
   //if(k>=800) cout<<"done searching."<<endl;                        
   for(int j=0;j<k;j++){
     int mynn_idx = nn_idx[j];
     //      Galaxy *g1 = galaxies[i];                                     
     //Galaxy *g2 = galaxies[mynn_idx];                                    
     if(fabs(galaxies[i]->Z()-galaxies[mynn_idx]->Z())<delta_z){
       ++nn;
     }
     if(nn==k_want){
       dist = galaxies[i]->ComovDist(sqrt(dists[j]));
     if(j>jmax) jmax=j;
     out<<dist<<" "<<j<<" "<<k<< endl;
     distfound = true;
     break;
   }
 }
 if(!distfound){
   dist = galaxies[i]->ComovDist(sqrt(dists[k-1]));
   if((dist<15)&&(k<n_pts/4.)){
     //cout<<dist<<" "<<i<<" ** "<<k<<endl;                              
     //re-search if you didn't look far out enough                       
     //only do this if your distance was smaller than the max of the data, and                                                                        
       //only do this if it won't put you over the max number of data points                                                                            
       k = k*2;
     nn_idx = new ANNidx[k];                     // allocate near neighindices                                                                       
  dists = new ANNdist[k];                     // allocate near neighbor dists                                                                      
  }
   else{
     out<<dist<<" 0 "<<k<< endl;
     distfound = true;
     break;
   }
 }
}
nndist[i] = dist;
}

delete [] nn_idx;
delete [] dists;
delete the_tree;

return nndist;
}
*/

#ifdef COLORS_FROM_RELATIVE_DENSITY
class GalaxyPercent{
public: 
  GalaxyPercent(float nndist, int gid):nndist(nndist), id(gid){};
  float dens()const{return nndist;};
  int gid()const{return id;};
  void dens(float tnndist){
    nndist = tnndist;
  }
  void gid(int tid){
    id = tid;
  }
private:
  float nndist;
  int id;
};

int CompareByDens (GalaxyPercent a, GalaxyPercent b)
{
  //return a->dens() > b->dens();
  return a.dens() > b.dens();
}

vector <float> GetNeighborPercents(vector <float> nndist, vector <Galaxy *> &galaxies)
{
  cout<<"normal version"<<endl;
  vector <float> color_percent(galaxies.size());
  cout<<"Made color percent vec"<<endl;
  int nColorBins = floor(zmax->GetVal()/ColorBinSize)+1;
  cout<<"got num of color bins"<<endl;
  struct GalaxyPercent tGalaxyId(0.,0);
  cout<<"made galaxy percent struct"<<endl;
  vector <GalaxyPercent> GalaxyId;
  cout<<"made galaxy id vec"<<endl;
  //  GalaxyId.reserve(galaxies.size());
  int nInColorBins = 0;
  int nInColorBin[nColorBins];
  cout<<"made n in color bin int vec"<<endl;
  for (int i=0;i<nColorBins;i++)
    nInColorBin[i] = 0;
  cout<<"set all to zero"<<endl;
  int min_gid = 1000000;
  int max_gid = 0;
  int min_i = 1000000;
  int max_i = 0;

  cout<<" Looping through percent bins."<<endl;
  for (int bin=0;bin<nColorBins;bin++){
    float binmin = bin*ColorBinSize;
    float binmax = binmin + ColorBinSize;
    //cout<<"  doing bin "<<bin<<" with min/max z = "<<binmin<<"/"<<binmax<<endl;
    for (int i=0;i<galaxies.size();i++){
      if (galaxies[i]->Z() < binmin || galaxies[i]->Z() > binmax)
	continue;
      float tnndist = nndist[i];
      tGalaxyId.dens(tnndist);
      tGalaxyId.gid(i);
      GalaxyId.push_back(tGalaxyId);
      nInColorBin[bin]++;
    }
    //cout<<"  sorting ID's."<<endl;
    /*
    cout<<"Before sort, first 3: "<<GalaxyId[0].gid()<<" "<<
      GalaxyId[0].dens()<<", "<<GalaxyId[1].gid()<<" "<<
      GalaxyId[1].dens()<<", "<<GalaxyId[2].gid()<<" "<<
      GalaxyId[2].dens()<<endl;
    */

    //check how many 0's are in the dens list
    int nZero = 0;
    for (int i=0;i<nInColorBin[bin];i++)
      if (GalaxyId[i].dens() <= 0.)
	nZero++;
    if (nZero > 0)
      cout<<" Before sorting there are "<<nZero<<" 0's in the anticipated range."<<endl;

    sort(GalaxyId.begin(),GalaxyId.end(),CompareByDens);

    //check how many 0's are in the dens list
    nZero = 0;
    for (int i=0;i<nInColorBin[bin];i++)
      if (GalaxyId[i].dens() <= 0.)
	nZero++;
    if (nZero > 0)
      cout<<" After sorting there are "<<nZero<<" 0's in the anticipated range."<<endl;

    /*
    cout<<"Sorted by dens: "<<GalaxyId[0].gid()<<" "<<
      GalaxyId[0].dens()<<", "<<GalaxyId[1].gid()<<" "<<
      GalaxyId[1].dens()<<", "<<GalaxyId[2].gid()<<" "<<
      GalaxyId[2].dens()<<endl;
    */

    //cout<<"  doing final percent calculation."<<endl;
    for (int i=0;i<nInColorBin[bin];i++){
      color_percent[GalaxyId[i].gid()] = float(i)/float(nInColorBin[bin]);
      if (GalaxyId[i].gid() > max_gid)
	max_gid = GalaxyId[i].gid();
      if (GalaxyId[i].gid() < min_gid)
	min_gid = GalaxyId[i].gid();
    }

    GalaxyId.clear();

    cout<<"     This bin had "<<nInColorBin[bin]<<" particles in it."<<endl;
    nInColorBins += nInColorBin[bin];
    //cout<<"  done with the bin."<<endl;
  }
  cout<<"First 3 color percents: "<<color_percent[0]<<" "<<color_percent[1]<<" "<<color_percent[2]<<endl;
  cout<<"Last 3 color percents: "<<color_percent[galaxies.size()-1]<<" "<<color_percent[galaxies.size()-2]<<" "<<color_percent[galaxies.size()-1]<<endl;
  cout<<"Did the percent for a total of "<<nInColorBins<<" Galaxies (tot # = "<<galaxies.size()<<")"<<endl;
  cout<<"Range of galaxy id's covered: "<<min_gid<<" - "<<max_gid<<endl;
  return color_percent;
}
#endif

// new 2d version
// this calculates the proper ra-dec distance
// Galaxies:  The n nearest objects to the points
// Points:  The locations from which the distances are measured
/*
vector <float> GetNeighborDist(vector <Galaxy *> galaxies, vector <Galaxy *> points){
  cout<<"In GetNeighborDist"<<endl;
  cout<<"  number of galaxies: "<<galaxies.size()<<endl;;
  cout<<"  number of points:   "<<points.size()<<endl;
  PRNTV(zmax->GetVal());
  //PRNTV(galaxies.size());
  int nbins = (int) (ceil(zmax->GetVal()/zbsize));
  PRNTV(nbins);
  vector <int> ginbin(nbins+1);
  for(int i=0;i<ginbin.size();i++){
    ginbin[i] = 0;
  }
  int OutNext = 0;
  cout<<" initial ginbin's: "<<ginbin[0]<<" "<<ginbin[1]<<" "<<ginbin[2]<<" "<<ginbin[3]<<" "<<ginbin[4]<<endl;
  for(int i=0;i<galaxies.size();i++){
    if(i>OutNext){
      //cout<<i<<" of "<<galaxies.size()<<endl;;
      OutNext += 1000;
    }
    ginbin[galaxies[i]->Zbin()]++;
    if (i < 5){
      cout<<galaxies[i]->Ra()<<" "<<galaxies[i]->Dec()<<" "<<galaxies[i]->Z()<<" "<<galaxies[i]->Zbin()<<" "<<ginbin[galaxies[i]->Zbin()]<<endl;
    }
  
  int ginanybin = 0;
  for(int i=0;i<ginbin.size();i++){
    ginanybin += ginbin[i];
    //cout<<"  "<<i<<" "<<ginbin[i]<<endl;
  }
  cout<<" Number of galaxies in bins (should only be ngals if Magmin >= Magmin_dens): "<<ginanybin<<endl;

  vector <int> pinbin(nbins+1);
  for(int i=0;i<pinbin.size();i++){
    pinbin[i] = 0;
  }
  cout<<"  Zeroed pinbin."<<endl;
  OutNext = 0;
  for(int i=0;i<points.size();i++){
    //cout<<i<<" of "<<galaxies.size()<<endl;
    if(i>OutNext){
      //cout<<i<<" of "<<galaxies.size()<<endl;;
      OutNext += 1000;
    }
    pinbin[points[i]->Zbin()]++;
    if (i < 5){
      cout<<points[i]->Ra()<<" "<<points[i]->Dec()<<" "<<points[i]->Z()<<" "<<points[i]->Zbin()<<" "<<pinbin[points[i]->Zbin()]<<endl;
    }
  }
  cout<<"  Set pinbin"<<endl;
  int pinanybin = 0;
  for(int i=0;i<pinbin.size();i++)
    pinanybin += pinbin[i];
  cout<<" Number of points in bins: "<<pinanybin<<endl;

#ifdef DEBUG
  cout<<"[neighbor] "<<nbins<<" z bins"<<endl;
#endif
  int startbin=0;
  bool have_started = false;
  for(int i=0;i<nbins;i++){
#ifdef DEBUG
        cout<<"[neighbor] "<<i<<" "<<ginbin[i]<<endl;
#endif
	//	if((ginbin[i]<10)&&(i>1)){
	if(ginbin[i]<10){
	  if(have_started == false){
	    startbin++;
	    //cout<<"not enough objects in bin "<<i<<endl;
	  }
	  else{
	    //cout<<"not enough objects in bin "<<i<<endl;
	    nbins = i;
	    break;
	  }
	}
	else have_started = true;
  }
  cout<<"using "<<startbin<<" "<<nbins<<endl;
  //vector <float> nndist(galaxies.size());
  vector <float> nndist(points.size());
  ANNpoint            query_pt;               // query point
  ANNidxArray         nn_idx;                 // near neighbor indices
  ANNdistArray        dists;                  // near neighbor distances
  ANNkd_tree * the_tree;
  double eps=0;
  int dim = 2;
  int k_want = 10;  // 10th nearest neighbor (first one is self)
  if(dens_measure == FIFTH)
    k_want = 5;

  query_pt = annAllocPt(dim);             // allocate query point
  //  int jmax = 0;
  ofstream ddout("dd2.out");
  ofstream nfout("nf.out");
  //loop over the z slices
  int  pi=0;
#ifdef DEBUG
  cout<<"Looping over "<<nbins-startbin<<"z bins"<<endl;
#endif
  int max_search= 350;

  for(int bi=startbin;bi<nbins;bi++){
    ANNpointArray       data_pts;        // data points to look for neighbors
    ANNpointArray       data_pts1;       // data points to measure dists for
    //int npts=ginbin[bi];
    int npts=pinbin[bi];
    if (npts == 0)
      continue;
    //cout<<"npts:"<<npts<<endl;
    int nsearch=ginbin[bi];
    //cout<<" bin has "<<npts<<" bright gals and "<<nsearch<<" total gals."<<endl;
    int k_search = 70;                   //initial search radius
    //if(npts<k_search) k_search=npts-1;
    if(nsearch<k_search) k_search=nsearch-1;
    nn_idx = new ANNidx[k_search];          // allocate near neigh indices
    dists = new ANNdist[k_search];          // allocate near neighbor dists

    data_pts1 = annAllocPts(npts, dim);         // allocate data points
    vector <int> ids(npts);
    vector <int> sids;
    sids.reserve(nsearch);
    if(bi>0) nsearch+=ginbin[bi-1];
    if(bi<ginbin.size()) nsearch+=ginbin[bi+1];
    data_pts = annAllocPts(nsearch, dim);         // allocate data points
    int nsi = 0;
    int npi = 0;

    //cout<<"Searching through bin "<<bi<<endl;
    //cout<<"  allocated nsearch = "<<nsearch<<endl;

//    for(int i=0;i<galaxies.size();i++){
//      if((galaxies[i]->Zbin()==bi)||
//	 (galaxies[i]->Zbin()==bi-1)||
//	 (galaxies[i]->Zbin()==bi+1)){
//	assert(nsi<nsearch);
//	data_pts[nsi][0] = galaxies[i]->Ra();
//	data_pts[nsi][1] = galaxies[i]->Dec();
//	sids.push_back(i);
//	nsi++;
//      }
//      if(galaxies[i]->Zbin()==bi){
//	data_pts1[npi][0] = galaxies[i]->Ra();
//	data_pts1[npi][1] = galaxies[i]->Dec();
//	ids[npi]=i;
//	npi++;
//      } 
//    }
    for(int i=0;i<galaxies.size();i++){
      if((galaxies[i]->Zbin()==bi)||
	 (galaxies[i]->Zbin()==bi-1)||
	 (galaxies[i]->Zbin()==bi+1)){
	assert(nsi<nsearch);
	data_pts[nsi][0] = galaxies[i]->Ra();
	data_pts[nsi][1] = galaxies[i]->Dec();
	sids.push_back(i);
	nsi++;
      }
    }
    for(int i=0;i<points.size();i++){
      if(points[i]->Zbin()==bi){
	data_pts1[npi][0] = points[i]->Ra();
	data_pts1[npi][1] = points[i]->Dec();
	ids[npi]=i;
	npi++;
      } 
    }

    //cout<<"The Tree."<<endl;
    the_tree = new ANNkd_tree(
			      data_pts,  // the data points
			      nsearch,   // number of points
			      dim);		      // dimension of space
    //cout<<"bi = "<<bi<<" of "<<nbins<<" npts = "<<npts<<" npi = "<<npi<<" nsi = "<<nsi<<endl;
    for(int i=0;i<npts;i++){
      //cout<<"Searching for "<<i<<" of "<<npts<<endl;
      int gid = ids[i];
      #ifdef DEBUG
      if(pi%50000==0) {cout<<bi<<" "<<i<<" "<<pi<<endl; system("date");}
      #endif
      bool distfound = false;
      ANNpoint query_pt = data_pts1[i];
      float dist=0;
      int n_found = -1; 
      int t_found = -1; 
      while(!distfound){
	//if(npts<k_search) k_search=npts-1;
	//if(nsi<k_search) k_search=npts-1;
	if(nsearch<k_search) k_search=nsearch-1;
	if(k_search>max_search) k_search=max_search;
	int nn = 0;
	//cout<<"doing search with k_search = "<<k_search<<endl;
	the_tree->annkSearch(			// search
			     query_pt,		// query point
			     k_search,		// number of near neighbors
			     nn_idx,		// nearest neighbors (returned)
			     dists,		// distance (returned)
			     eps);		// error bound  
	for(int j=0;j<k_search;j++){
	  int chosen_id = sids[nn_idx[j]];
	  //  if((gid!=chosen_id)
	  //&&(fabs(galaxies[gid]->Z()-galaxies[chosen_id]->Z())<delta_z)){
	  //if((fabs(galaxies[gid]->Z()-galaxies[chosen_id]->Z())<delta_z)){
	  if((fabs(points[gid]->Z()-galaxies[chosen_id]->Z())<delta_z)){
	    ++nn;
	  }
	  if(nn==k_want){	
	    //for(int jjj=0;jjj<j;jjj++){
	    //cout<<galaxies[gid]->ComovDist(sqrt(dists[jjj]))<<" ";
	    //}
	    //cout<<endl;
	    //dist = galaxies[gid]->ComovDist(sqrt(dists[j]));
	    dist = points[gid]->ComovDist(sqrt(dists[j]));
	    n_found = j;
	    t_found = 0;
	    distfound = true;
	    break;
	  }
	}
	if(!distfound){
	  //dist = galaxies[pi]->ComovDist(sqrt(dists[k_search-1]));
	  dist = points[pi]->ComovDist(sqrt(dists[k_search-1]));
	  //if you're at the end, pretend you've found it
	  //if((k_search==npts-1)||(k_search>=max_search)){
	  if((k_search==nsearch-1)||(k_search>=max_search)){
	    n_found = k_search-1;
	    t_found = 1;
	    distfound = true;
	  }
	  //if you're big already, pretend you've found it
	  //	  if(dist>15){
	  if(dist>10){
	    n_found = -1;
	    t_found = 2;
	    distfound = true;
	  }
	  else{
	    //nfout<<dist<<" "<<i<<" "<<galaxies[i]->Z()<<" "<<k_search<<endl;
	    nfout<<dist<<" "<<i<<" "<<points[i]->Z()<<" "<<k_search<<endl;
	    //re-search if you didn't look far out enough
	    k_search = 2*k_search;
	    delete [] nn_idx;
	    delete [] dists;
	    nn_idx = new ANNidx[k_search];
	    dists = new ANNdist[k_search];
	  }
  	}//if dist is not found
      }//while dist is not found
      nndist[gid] = dist;
      ddout<<i<<" "<<pi<<" "<<gid<<" "<<dist<<" "<<t_found<<" "<<n_found<<" "<<k_search<<" "<<npts<<" "<<bi<<endl;
      pi++;
    }//loop over points in z bin
    //    ofstream ddout("dd2.out");
    //nnout<<nndist[i]<<endl;
    delete the_tree;
    delete [] nn_idx;
    delete [] dists;  
  }//loop over bins
  //MSG("Exiting neighbor");
  return nndist;
}
*/

/*
// slow 2d version
// this calculates the proper ra-dec distance
vector <float> GetNeighborDistOrig(vector <Galaxy *> galaxies){
  vector <float> nndist(galaxies.size());
  int                 n_pts = galaxies.size();                  // actual number of data points
  ANNpointArray       data_pts;               // data points
  ANNpoint            query_pt;               // query point
  ANNidxArray         nn_idx;                 // near neighbor indices
  ANNdistArray        dists;                  // near neighbor distances
  ANNkd_tree          *the_tree;              // search structure
  double eps=0;
  int dim = 2;
  int k_want = 11;  //counts itself, so this is really the 10th nearest neighbor
  if(dens_measure == FIFTH)
    k_want = 5;
  int k = 100;    
  query_pt = annAllocPt(dim);                 // allocate query point
  data_pts = annAllocPts(n_pts, dim);         // allocate data points
  nn_idx = new ANNidx[k];                     // allocate near neigh indices
  dists = new ANNdist[k];                     // allocate near neighbor dists
  
  int jmax = 0;

  //ofstream outpts("pts.out");
  for(int i=0;i<n_pts;i++){
    Particle *p = galaxies[i]->P();
    data_pts[i][0] = p->Ra();
    data_pts[i][1] = p->Dec();
    //data_pts[i][0] = p->X();
    //data_pts[i][1] = p->Y();
  //outpts<<data_pts[i][0]<<" "<<data_pts[i][1]<<endl;
  }

  ofstream outall("all.out");
  the_tree = new ANNkd_tree(			// build search structure
			    data_pts,			// the data points
			    n_pts,			// number of points
			    dim);			// dimension of space

  ofstream out("dd2.out");

  for(int i=0;i<n_pts;i++){
    if(i%20000==0) {cout<<i<<endl; system("date");}
    bool distfound = false;
 
    ANNpoint query_pt = data_pts[i];
    float dist=0;
    k=100;
    while(!distfound){
      int nn = 0;
      //if(k>=800) cout<<"searching..."<<k<<" ";
      the_tree->annkSearch(			// search
			   query_pt,		// query point
			   k,			// number of near neighbors
			   nn_idx,		// nearest neighbors (returned)
			   dists,		// distance (returned)
			   eps);		// error bound  
      //if(k>=800) cout<<"done searching."<<endl; 
      for(int j=0;j<k;j++){
	int mynn_idx = nn_idx[j];
	//	Galaxy *g1 = galaxies[i];
	//Galaxy *g2 = galaxies[mynn_idx];
	if(fabs(galaxies[i]->Z()-galaxies[mynn_idx]->Z())<delta_z){
	  ++nn;
	}
	if(nn==k_want){	
	  dist = galaxies[i]->ComovDist(sqrt(dists[j]));
	  //	cout<<dist<<" "<<galaxies[i]->ComovDist(dist)<<" "<<galaxies[i]->P()->Dist2D(galaxies[mynn_idx]->P())<<endl;
 	  if(j>jmax) jmax=j;
	  out<<dist<<" "<<j<<" "<<k<< endl;
	  //if(k!=100)cout<<dist<<" "<<i<<" "<<j<<" "<<k<<endl;
	  //if(dist==0)
	  //cout<<query_pt[0]<<" "<<data_pts[mynn_idx][0]<<" "<<j<<" "<<nn<<endl;
	  distfound = true;
	  break;
	}
      }
      if(!distfound){
	dist = galaxies[i]->ComovDist(sqrt(dists[k-1]));
	if((dist<15)&&(k<n_pts/4.)){
	  //cout<<dist<<" "<<i<<" ** "<<k<<endl;
	  //re-search if you didn't look far out enough
	  //only do this if your distance was smaller than the max of the data, and 
	  //only do this if it won't put you over the max number of data points
	  k = k*2;
	  delete [] nn_idx;
	  delete [] dists;
	  nn_idx = new ANNidx[k];                     // allocate near neigh indices
	  dists = new ANNdist[k];                     // allocate near neighbor dists
	}
	else{
	  out<<dist<<" 0 "<<k<< endl;
	  distfound = true;
	  break;
	}	  
      }
    }
    nndist[i] = dist;
  }

  delete [] nn_idx;
  delete [] dists;
  delete the_tree;

  return nndist;
}
*/


//vector <int> GetSEDs(vector <Galaxy *> &galaxies, vector <float> &nndist, vector <GalSED> & galseds){
vector <int> GetSEDs(vector <Galaxy *> &galaxies, vector <float> &nndist, vector <GalSED> & galseds, vector <Halo *> &halos){
  vector <int> sed_ids(galaxies.size());
  float Percent = 0.0;
  cout<<galaxies[0]->Mr()<<endl;
  //cout<<"Need to check some magnitudes and densities: "<<endl;
  //for(int i=-22;i<=-4;i++){
  //cout<<"  "<<i<<", "<<LumNumberDensityInterp(float(i))<<endl;
  //}
  for(int gi=0;gi<galaxies.size();gi++){
    //cout<<gi<<" "<<galaxies[gi]->Mr()<<endl;
    if(float(gi)/float(galaxies.size()) > Percent){
      cout<<"  "<<Percent*100.<<"% done"<<endl;
      Percent += 0.1;
    }
    //cout<<gi<<" of "<<galaxies.size()<<endl;

    //if(galaxies[gi]->Central())
    //{
    //extern double normal_random(float mean, float stddev);
    //Galaxy * gal = galaxies[gi];
    //Particle * p = gal->P();
    //assert(p);  //this better be true since you removed the other ones.
    //int hid = p->Hid();
    //double m200 = halos[hid]->M();
    ////cout<<"lookin at central galaxy ("<<galaxies[gi]->Central()<<") with m200="<<m200<<", density = "<<nndist[gi]<<", and original Mr = "<<galaxies[gi]->Mr();
    //double m14 = m200/1e14;
    //double lum = 4*pow(m14,0.3);
    //lum = pow(10.0,normal_random(log10(lum),0.15));
    //double mr = -2.5*log10(lum) - 20.44;
    //galaxies[gi]->Mr(mr);
    //nndist[gi] *= 0.5; //We want to make the BCGs dense
    ////cout<<" adjusted to "<<galaxies[gi]->Mr()<<endl;
    ////if(mr>Magmin_col) mr = Magmin_col;
    //}
    Particle * p = galaxies[gi]->P();
    assert(p);  //this better be true since you removed the other ones.
    //int hid = p->Hid();
    //if(hid>=0){
#ifdef DEBUG
    if(gi%20000==0) {cout<<gi<<endl; system("date");}
#endif
    double mr = galaxies[gi]->Mr();
    //if(mr>Magmin_col) mr = Magmin_col;
    //cout<<"non central galaxy with Mr = "<<galaxies[gi]->Mr()<<", dens = "<<nndist[gi]<<endl;
    //cout<<"Choosnig..."<<endl;

	sed_ids[gi] = findCloseGalaxies2(galseds, mr, nndist[gi], galaxies[gi]->Z(), galaxies[gi]->Central());
    //cout<<"Simulated galaxy: mr = "<<mr<<" nndist = "<<nndist[gi]<<endl;
    //cout<<"SED galaxy: mr = "<<galseds[sed_ids[gi]].MR()<<" nndist = "<<galseds[sed_ids[gi]].Dens()<<endl;
    //cout<<sed_ids[gi]<<endl;



	//sed_ids[gi] = ChooseSED(galseds,mr,nndist[gi], galaxies[gi]->Z(), galaxies[gi]->Central());
    
	
	
	//}
    //else sed_ids[gi]=-1;
  }
  return sed_ids;
}

//This is for use without neighbor distances
vector <int> GetSEDs(vector <Galaxy *> &galaxies, vector <GalSED> & galseds){
  vector <int> sed_ids(galaxies.size());
  for(int gi=0;gi<galaxies.size();gi++){
    Particle * p = galaxies[gi]->P();
    assert(p);  //this better be true since you removed the other ones.
    //int hid = p->Hid();
    //if(hid>=0){
#ifdef DEBUG
    if(gi%20000==0) {cout<<gi<<endl; system("date");}
#endif
    double mr = galaxies[gi]->Mr();
    //if(mr>Magmin_col) mr = Magmin_col;

    sed_ids[gi] = ChooseSED(galseds,mr,0, galaxies[gi]->Z(), galaxies[gi]->Central());
    //}
    //else sed_ids[gi]=-1;
  }
  return sed_ids;
}

//This is for use without neighbor distances
vector <int> GetBCG_SEDs(vector <Galaxy *> &galaxies, vector <GalSED> & galseds){
  vector <int> sed_ids(galaxies.size());
  for(int gi=0;gi<galaxies.size();gi++){

      //}
    //else sed_ids[gi]=-1;
  }
  return sed_ids;
}

