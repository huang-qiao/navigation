#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cstdio>

#include "map.hpp"

// Create a new map
MapPtr CreateMap(void)
{
  MapPtr map = std::make_shared<Map>();

  //map = (Map*) malloc(sizeof(Map));

  // Assume we start at (0, 0)
  map->origin_x = 0;
  map->origin_y = 0;

  // Make the size odd
  map->size_x = 0;
  map->size_y = 0;
  map->scale = 0;

  // Allocate storage for main map
  //map->cells = (MapCell*) NULL;
  map->cells.clear();

  return map;
}

// Destroy a map
void DeleteMap(MapPtr map)
{
  //free(map->cells);
  //free(map);
  map.reset();
  return;
}

// Get the cell at the given point
CellPtr Map::getCell(double ox, double oy, double oa)
{
  int i, j;
  CellPtr cell;

  i = toGridX(ox);
  j = toGridY(oy);

  if (!isValid(i, j))
    return NULL;

  //cell = map->cells + MAP_INDEX(map, i, j);
  cell = cells[toIndex(i, j)];

  return cell;
}
