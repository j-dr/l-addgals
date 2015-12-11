#ifndef SHAPES_H
#define SHAPES_H 1

#include <vector>
#include <iterator>
#include <istream>
#include <utility>
#include <functional>
#include "shapeconfig.h"

//extern prefstruct prefs;

struct shapemag 
{
  double bands[5];
};

std::istream & operator>>(std::istream & is, shapemag & in);

void generate_shapes(std::vector<double>& mags, std::vector<double>& e, std::vector<double>& s,
		     int nelem, int vl);

#endif
