#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "map.hpp"

using namespace amcl::map;

// Create a new map
//map_t *map_alloc(void)
MapPtr Map::CreateMap()
{
  MapPtr map = std::make_shared<Map>();

  // Assume we start at (0, 0)
  map->origin_x = 0;
  map->origin_y = 0;
  
  // Make the size odd
  map->size_x = 0;
  map->size_y = 0;
  map->scale = 0;
    
  return map;
}


/* Destroy a map
void map_free(map_t *map)
{
  free(map->cells);
  free(map);
  return;
}
*/

// Get the cell at the given point
//map_cell_t *map_get_cell(map_t *map, double ox, double oy, double oa)
CellPtr Map::getCell(const double &ox, const double &oy, const double &oa)
{
  int i = (int)toGridX(ox);
  int j = (int)toGridY(oy);

  if (!isValid(i, j))
    return nullptr;

  CellPtr cell = cells[toIndex(i, j)];
  return cell;
}

