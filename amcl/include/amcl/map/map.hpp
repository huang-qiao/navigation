#pragma once

#include <memory>
#include <vector>

// Limits
#define MAP_WIFI_MAX_LEVELS 8


// Description for a single map cell.
typedef struct
{
  // Occupancy state (-1 = free, 0 = unknown, +1 = occ)
  int occ_state;

  // Distance to the nearest occupied cell
  double occ_dist;

  // Wifi levels
  //int wifi_levels[MAP_WIFI_MAX_LEVELS];

} MapCell;

using Cell = MapCell;
using CellPtr = std::shared_ptr<Cell>;

// Description for a map
struct Map : std::enable_shared_from_this<Map>
{
  // Map origin; the map is a viewport onto a conceptual larger map.
  double origin_x, origin_y;

  // Map scale (m/cell)
  double scale;

  // Map dimensions (number of cells)
  int size_x, size_y;

  // The map data, stored as a grid
  //CellPtr cells;
  std::vector<CellPtr> cells;

  // Max distance at which we care about obstacles, for constructing
  // likelihood field
  double max_occ_dist;

  /**************************************************************************
   * Basic map functions
   **************************************************************************/

  // Get the cell at the given point
  CellPtr getCell(double ox, double oy, double oa);

  // Load an occupancy map
  int loadOcc(const char *filename, double scale, int negate);

  // Update the cspace distances
  void updateCSpace(double max_occ_dist);

  /**************************************************************************
   * Range functions
   **************************************************************************/

  // Extract a single range reading from the map
  double calcRange(double ox, double oy, double oa, double max_range);

  /**************************************************************************
   * Map manipulation
   **************************************************************************/

  // Convert from map index to world coords
  inline double toWorldX(double i) { return origin_x + ((i) - size_x / 2) * scale; }
  inline double toWorldY(double j) { return origin_y + ((j) - size_y / 2) * scale; }

  // Convert from world coords to map coords
  inline double toGridX(double x) { return floor((x - origin_x) / scale + 0.5) + size_x / 2; }
  inline double toGridY(double y) { return floor((y - origin_y) / scale + 0.5) + size_y / 2; }

  // Test to see if the given map coords lie within the absolute map bounds.
  inline bool isValid(int i, int j) { return (i >= 0) && (i < size_x) && (j >= 0) && (j < size_y); }

  // Compute the cell index for the given map coords.
  inline int toIndex(int i, int j) { return (i) + (j) * size_x; }
};

using MapPtr = std::shared_ptr<Map>;

/**************************************************************************
 * Global map functions
 **************************************************************************/

// Create a new (empty) map
MapPtr CreateMap(void);

// Destroy a map
void DeleteMap(MapPtr map);



// Load a wifi signal strength map
//int map_load_wifi(MapPtr map, const char *filename, int index);

// [TODO] re-enable map-drawing functions
//        replace rtk stuff with opencv
/**************************************************************************
 * GUI/diagnostic functions
 **************************************************************************/

// Draw the occupancy grid
void map_draw_occ(MapPtr map, struct _rtk_fig_t *fig);

// Draw the cspace map
void map_draw_cspace(MapPtr map, struct _rtk_fig_t *fig);

// Draw a wifi map
void map_draw_wifi(MapPtr map, struct _rtk_fig_t *fig, int index);

