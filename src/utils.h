#ifndef utils_h
#define utils_h
#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

struct data_point {
  float x, y, err;
};

using namespace std;

void read_data(string filename, data_point data[], int npoints);
void read_cdata(string filename, vector < vector <float> > &data, int npoints);
void	spline(float x[], float y[] ,int n, float yp1, float ypn, float y2[]);
void	splint(float xa[], float ya[], float y2a[], int n, float x, float *y);
void	spldint(float xa[], float ya[], float y2a[], int n, float x, float *y);
double	gammafn2(float a, float xmin);
double	gammln(float xx);
double	gammafn1(float z);

#endif
