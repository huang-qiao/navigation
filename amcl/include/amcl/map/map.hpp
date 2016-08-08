#pragma once

#include <cmath>
#include <cstdint>
#include <memory>
#include <vector>
#include <queue>

struct CellData;
using CellDataPtr = std::shared_ptr<CellData>;

bool operator<(const CellData& a, const CellData& b);

class CachedDistanceMap;
using CachedDistanceMapPtr = std::shared_ptr<CachedDistanceMap>;


namespace amcl {
namespace map {

// forward declaration
struct Cell;
using CellPtr = std::shared_ptr<Cell>;


struct Map;
using MapPtr = std::shared_ptr<Map>;

// Description for a single map cell.
struct Cell
{
  // Occupancy state (-1 = free, 0 = unknown, +1 = occ)
  int occ_state;
  // Distance to the nearest occupied cell
  double occ_dist;
};

// Description for a map
struct Map
{
  // Create a new map
  static MapPtr CreateMap();

  // Load an occupancy grid
  static MapPtr LoadOccMap(const std::string &filename, const double &scale, const int &negate);

  // Map origin; the map is a viewport onto a conceptual larger map.
  double origin_x, origin_y;
  
  // Map scale (m/cell)
  double scale;

  // Map dimensions (number of cells)
  int size_x, size_y;
  
  // The map data, stored as a grid
  std::vector<CellPtr> cells;

  // Max distance at which we care about obstacles, for constructing
  // likelihood field
  double max_occ_dist_;

  CellPtr getCell(const double &ox, const double &oy, const double &oa = 0.0);

  void updateCSpace(const double &max_occ_dist);

  // Extract a single range reading from the map.  Unknown cells and/or
  // out-of-bound cells are treated as occupied, which makes it easy to
  // use Stage bitmap files.
  double calcRange(const double &ox, const double &oy, const double &oa, const double &max_range);

  inline double toWorldX(const double &gx) { return origin_x + ((gx) - size_x / 2) * scale; }
  inline double toWorldY(const double &gy) { return origin_y + ((gy) - size_y / 2) * scale; }
  inline double toGridX(const double &wx) { return floor((wx - origin_x) / scale + 0.5) + size_x / 2; }
  inline double toGridY(const double &wy) { return floor((wy - origin_y) / scale + 0.5) + size_y / 2; }

  // Test to see if the given map coords lie within the absolute map bounds.
  inline bool isValid(double i, double j) { return (i >= 0) && (i < size_x) && (j >= 0) && (j < size_y); }
  // Compute the cell index for the given map coords.
  inline size_t toIndex(size_t i, size_t j) { return i+j*size_x; }

private:
  void enqueue(const int &i, const int &j, const int &src_i, const int &src_j, std::priority_queue<CellData> &Q, CachedDistanceMapPtr cdm, std::vector<bool> &marked);
};


/**************************************************************************
 * Basic map functions
 **************************************************************************/

// Create a new (empty) map
//map_t *map_alloc(void);

// Destroy a map
//void map_free(map_t *map);

// Get the cell at the given point
//map_cell_t* map_get_cell(map_t *map, double ox, double oy, double oa);

// Load an occupancy map
//int map_load_occ(map_t *map, const char *filename, double scale, int negate);

// Load a wifi signal strength map
//int map_load_wifi(map_t *map, const char *filename, int index);

// Update the cspace distances
//void map_update_cspace(map_t *map, double max_occ_dist);


/**************************************************************************
 * Range functions
 **************************************************************************/

// Extract a single range reading from the map
//double map_calc_range(map_t *map, double ox, double oy, double oa, double max_range);


/**************************************************************************
 * GUI/diagnostic functions
 **************************************************************************/

// Draw the occupancy grid
//void map_draw_occ(map_t *map, struct _rtk_fig_t *fig);

// Draw the cspace map
//void map_draw_cspace(map_t *map, struct _rtk_fig_t *fig);

// Draw a wifi map
//void map_draw_wifi(map_t *map, struct _rtk_fig_t *fig, int index);


/**************************************************************************
 * Map manipulation macros
 **************************************************************************

// Convert from map index to world coords
#define MAP_WXGX(map, i) (map->origin_x + ((i) - map->size_x / 2) * map->scale)
#define MAP_WYGY(map, j) (map->origin_y + ((j) - map->size_y / 2) * map->scale)

// Convert from world coords to map coords
#define MAP_GXWX(map, x) (floor((x - map->origin_x) / map->scale + 0.5) + map->size_x / 2)
#define MAP_GYWY(map, y) (floor((y - map->origin_y) / map->scale + 0.5) + map->size_y / 2)

// Test to see if the given map coords lie within the absolute map bounds.
#define MAP_VALID(map, i, j) ((i >= 0) && (i < map->size_x) && (j >= 0) && (j < map->size_y))

// Compute the cell index for the given map coords.
#define MAP_INDEX(map, i, j) ((i) + (j) * map->size_x)
*/
} // namespace map
} // namespace amcl
