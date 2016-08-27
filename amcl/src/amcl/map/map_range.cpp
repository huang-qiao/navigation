#include <cassert>
#include <cmath>
#include <cstring>
#include <cstdlib>

#include "map.hpp"

// Extract a single range reading from the map.  Unknown cells and/or
// out-of-bound cells are treated as occupied, which makes it easy to
// use Stage bitmap files.
double Map::calcRange(double ox, double oy, double oa, double max_range)
{
  // Bresenham raytracing
  int x0,x1,y0,y1;
  int x,y;
  int xstep, ystep;
  char steep;
  int tmp;
  int deltax, deltay, error, deltaerr;

  x0 = toGridX(ox);
  y0 = toGridY(oy);

  x1 = toGridX(ox + max_range * cos(oa));
  y1 = toGridY(oy + max_range * sin(oa));

  if(abs(y1-y0) > abs(x1-x0))
    steep = 1;
  else
    steep = 0;

  if(steep)
  {
    tmp = x0;
    x0 = y0;
    y0 = tmp;

    tmp = x1;
    x1 = y1;
    y1 = tmp;
  }

  deltax = abs(x1-x0);
  deltay = abs(y1-y0);
  error = 0;
  deltaerr = deltay;

  x = x0;
  y = y0;

  if(x0 < x1)
    xstep = 1;
  else
    xstep = -1;
  if(y0 < y1)
    ystep = 1;
  else
    ystep = -1;

  if(steep)
  {
    if(!isValid(y,x) || cells[toIndex(y,x)]->occ_state > -1)
      return sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0)) * scale;
  }
  else
  {
    if(!isValid(x,y) || cells[toIndex(x,y)]->occ_state > -1)
      return sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0)) * scale;
  }

  while(x != (x1 + xstep * 1))
  {
    x += xstep;
    error += deltaerr;
    if(2*error >= deltax)
    {
      y += ystep;
      error -= deltax;
    }

    if(steep)
    {
      if(!isValid(y,x) || cells[toIndex(y,x)]->occ_state > -1)
        return sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0)) * scale;
    }
    else
    {
      if(!isValid(x,y) || cells[toIndex(x,y)]->occ_state > -1)
        return sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0)) * scale;
    }
  }
  return max_range;
}
