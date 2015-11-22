#ifndef KDTREE_H
#define KDTREE_H
#include <vector>
#include "particle.h"

template <int N, typename T> class KDNode
{
 public:
  KDNode();
  KDNode(T value);
  float point[N];
  KDNode *left;
  KDNode *right;
}

template <int N, typename T> bool compareNodes(KDNode<N,T> lhs, KDNode<N,T> rhs, int dim);

template <int N, typename T> class KDTree
{
 public:
  KDTree();
  KDTree(std::vector<T> points);
  void insert(T pt);
  KDNode<N, T>& findNode(const T search);
  T& nearestNeighbor(const T search);

  int numnodes;
  int depth;
  KDNode<N, T> *treenodes;
}

template <int N, typename T> KDNode<N,T>& makeTree(KDNode &t, int len, int dim);
template <int N, typename T> KDNode<N,T>& findMedian(KDNode<N,T> &start, KDNode<N,T> &end, int idx);

#endif

