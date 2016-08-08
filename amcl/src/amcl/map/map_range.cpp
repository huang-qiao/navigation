#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "map.hpp"

using namespace amcl::map;

// Extract a single range reading from the map.  Unknown cells and/or
// out-of-bound cells are treated as occupied, which makes it easy to
// use Stage bitmap files.
//double map_calc_range(map_t *map, double ox, double oy, double oa, double max_range)
double Map::calcRange(const double &ox, const double &oy, const double &oa, const double &max_range)
{
  // Bresenham raytracing

  // start cell
  int x0 = (int)toGridX(ox);
  int y0 = (int)toGridY(oy);
  
  // end cell
  int x1 = (int)toGridX(ox + max_range * cos(oa));
  int y1 = (int)toGridY(oy + max_range * sin(oa));

  bool steep = (abs(y1-y0) > abs(x1-x0))? true : false;

  // if steep, switch dX, dY
  // [FIXME] it's kinda ugly, is it necessary?
  if(steep) {
    int tmp;
    tmp = x0;
    x0 = y0;
    y0 = tmp;

    tmp = x1;
    x1 = y1;
    y1 = tmp;
  }

  int deltax = abs(x1-x0);
  int deltay = abs(y1-y0);
  int error = 0;
  int deltaerr = deltay;

  // [FIXME] very confusing here. since x=x0, y=y0, then
  //         it always return 0 if cell is not valid or
  //         cell=[unknown, occupied]
  int x = x0;
  int y = y0;

  int xstep = (x0 < x1)? 1 : -1;
  int ystep = (y0 < y1)? 1 : -1;


  if(steep) {
    // x,y are switched
    if(!isValid(y,x) || cells[toIndex(y,x)]->occ_state > -1)
      return sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0)) * scale;
  } else {
    if(!isValid(x,y) || cells[toIndex(x,y)]->occ_state > -1)
      return sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0)) * scale;
  }
  // [FIXME] very confusing here. since x=x0, y=y0, then
  //         it always return 0 if cell is not valid or
  //         cell=[unknown, occupied]

  while(x != (x1 + xstep * 1)) {
    x += xstep;
    error += deltaerr;
    if(2*error >= deltax) {
      y += ystep;
      error -= deltax;
    }

    if(steep) {
      if(!isValid(y,x) || cells[toIndex(y,x)]->occ_state > -1)
        return sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0)) * scale;
    } else {
      if(!isValid(x,y) || cells[toIndex(x,y)]->occ_state > -1)
        return sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0)) * scale;
    }
  }
  return max_range;
}
