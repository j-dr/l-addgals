#include "kdtree.h"
#include <vector>


template <int N, typename T> bool compareNodes(int dim, KDNode<N,T> &lhs, KDNode<N,T> &rhs)
{
  return lhs.point[dim] < rhs.point[dim];
}

inline template<int N, typename T> void swap(KDNode &x, KDNode &y) {
  float tmp[N];
  memcpy(tmp,  x->point, sizeof(tmp));
  memcpy(x->point, y->point, sizeof(tmp));
  memcpy(y->point, tmp,  sizeof(tmp));
}
  
template <int N, typename T> KDNode<N, T>& findMedian(KDNode<N,T>& start, KDNode<N,T> &end, int idx)
{
  if (end <= start) return NULL;
  if (end == start + 1)
    return start;
 
  KDNode<N,T> *p, *store, *md = start + (end - start) / 2;
  float pivot;

  while (1) {
    pivot = md->x[idx];
 
    swap(md, end - 1);
    for (store = p = start; p < end; p++) {
      if (p->x[idx] < pivot) {
	if (p != store)
	  swap(p, store);
	store++;
      }
    }
    swap(store, end - 1);
 
    /* median has duplicate values */
    if (store->x[idx] == md->x[idx])
      return md;
 
    if (store > md) end = store;
    else        start = store;
  }
}


  
template <int N, typename T> KDNode<N, T>& makeTree(KDNode<N,T>& t, int len, int dim)
{

  KDNode<N,T> *n;
  if (!len) return NULL;
 
  if ((n = find_median(t, t + len, dim))) {
    dim = (dim + 1) % N;
    n->left  = make_tree(t, n - t, dim);
    n->right = make_tree(n + 1, t + len - (n + 1), dim);
  }
  return n;
}

template <int N, typename T> KDTree::KDTree(vector<T> points)
{
  //Create Nodes from points
  KDNode<N,T> *treenodes;
  vector<KDNode<N,T>> nodes(points.size());
  vector<vector <T>> sortednodes(N);
  int i;

  for (vector<T>::iterator itr=points.begin(); itr!=points.end(); itr++)
    {
      nodes[(itr-points.begin())] = KDNode(*itr);
    }

  for (i=0; i<N; i++)
    {
      copy(nodes.begin(), nodes.end(), sortednodes[i].begin());
      sort(sortednodes[i].begin(), sortednodes[i].end()
	   bind(compareNodes<N,T>, i, _2, _3));
    }

  treenodes = makeTree<N, T> (&nodes[0], nodes.size(), 0);
}

  
