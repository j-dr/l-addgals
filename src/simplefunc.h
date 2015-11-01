#ifndef simplefunc_h
#define simplefunc_h

#include <cmath>
#include <iostream>
using namespace std;

template <class T> inline T sqr(T x) {return (x)*(x);}
template <class T> inline T cube(T x) { return (x)*(x)*(x);}
template <class T> inline T cuberoot(T x) { return pow((x),1/.3);}
template <class T> inline T abs(T x) {
  if (x>=0) return x;
  else return 0-x;
}

//#define max(a,b) ((a)>(b)?(a):(b))
//#define min(a,b) ((a)<(b)?(a):(b))
#define PRNTV(x) std::cout << #x " = " << x << "\n";
#define PRNT(f,x) std::cout <<"[" <<f<<"] "<< #x " = " << x << "\n";
#define MSG(x) std::cout <<x<<endl;
#define PRNTVS(f,x,s) std::cout <<"[" <<f<<"] "<< #x " = " << x <<s<< "\n";


#endif
