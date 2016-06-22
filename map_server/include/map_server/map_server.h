#pragma once

#include <string>
#include <fstream>
#include <memory>
#include <opencv2/opencv.hpp>
namespace map_server {
/** Map mode
 *  Default: TRINARY -
 *      value >= occ_th - Occupied (100)
 *      value <= free_th - Free (0)
 *      otherwise - Unknown
 *  SCALE -
 *      alpha < 1.0 - Unknown
 *      value >= occ_th - Occupied (100)
 *      value <= free_th - Free (0)
 *      otherwise - f( (free_th, occ_th) ) = (0, 100)
 *          (linearly map in between values to (0,100)
 *  RAW -
 *      value = value
 */
enum MapMode { TRINARY, SCALE, RAW };

struct OccupancyGridMap
{
  // meta data
  long map_load_time;
  double resolution; // [m/cell]
  unsigned int width; // [cells]
  unsigned int height; // [cells]
  double origin_x, origin_y, origin_th; // [m,m,rad]
  // map data
  std::vector<int> data;
};

using GridMap = OccupancyGridMap;

bool LoadMapFromFile(const std::string& fname, std::shared_ptr<GridMap> gridmap);

/** Read the image from file and fill out the resp object, for later
 * use when our services are requested.
 *
 * @param resp The map wil be written into here
 * @param fname The image file to read from
 * @param res The resolution of the map (gets stored in resp)
 * @param negate If true, then whiter pixels are occupied, and blacker
 *               pixels are free
 * @param occ_th Threshold above which pixels are occupied
 * @param free_th Threshold below which pixels are free
 * @param origin Triple specifying 2-D pose of lower-left corner of image
 * @param mode Map mode
 * @throws std::runtime_error If the image file can't be loaded
 **/
bool LoadMapFromFile(const char* fname, double res, bool negate, double occ_th, double free_th, double* origin, MapMode mode = TRINARY, std::shared_ptr<GridMap> map = nullptr);

} // namespace map_server
