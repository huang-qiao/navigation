#include <cstring>
#include <stdexcept>

#include <cstdlib>
#include <cstdio>
#include <string>
#include <memory>

#include <opencv2/opencv.hpp>
#include "map_server.h"

// compute linear index for given map coords
#define MAP_IDX(sx, i, j) ((sx) * (j) + (i))

#include <cstring>
#define __FILENAME__ \
  (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)

#define TRACE_FUNC \
  do { \
    std::cout << __FILE__ << ":" << __LINE__ << " in " << __func__ << std::endl; \
  } while (0);


namespace map_server {

bool
LoadMapFromFile(const std::string& fname, std::shared_ptr<GridMap> gridmap)
{
  std::string mapfname;
  double origin[3];
  int negate;
  double res, occ_th, free_th;
  MapMode mode = TRINARY;

  std::ifstream fin(fname.c_str());
  if (fin.fail()) {
    std::cout << __FILENAME__ << ":" << __LINE__ << " in " << __func__ << "could not open " << fname.c_str();
    exit(-1);
  }
  cv::FileStorage doc(std::string(fname), cv::FileStorage::READ);
  //cv::FileNode doc;
  try {
    doc["resolution"] >> res;
  } catch (std::exception e) {
    std::cout << e.what() << std::endl;
    std::cout << "The map does not contain a resolution tag or it is invalid.";
    exit(-1);
  }
  try {
    doc["negate"] >> negate;
  } catch (std::exception e) {
    std::cout << "The map does not contain a negate tag or it is invalid.";
    exit(-1);
  }
  try {
    doc["occupied_thresh"] >> occ_th;
  } catch (std::exception e) {
    std::cout << e.what();
    std::cout
        << "The map does not contain an occupied_thresh tag or it is invalid.";
    exit(-1);
  }
  try {
    doc["free_thresh"] >> free_th;
  } catch (std::exception e) {
    std::cout << e.what();
    std::cout << "The map does not contain a free_thresh tag or it is invalid.";
    exit(-1);
  }
  try {
    std::string modeS = "";
    doc["mode"] >> modeS;
    if (modeS == "trinary")
      mode = TRINARY;
    else if (modeS == "scale")
      mode = SCALE;
    else if (modeS == "raw")
      mode = RAW;
    else {
      std::cerr << "Invalid mode tag " << modeS.c_str() << std::endl;
      std::cerr << "Set TRINARY as defult" << std::endl;
      mode = TRINARY;
    }
  } catch (std::exception e) {
    std::cout << e.what();
    std::cout << "The map does not contain a mode tag or it is invalid... "
                 "assuming Trinary";
    mode = TRINARY;
  }
  try {
    doc["origin"][0] >> origin[0];
    doc["origin"][1] >> origin[1];
    doc["origin"][2] >> origin[2];
  } catch (std::exception e) {
    std::cout << e.what();
    std::cout << "The map does not contain an origin tag or it is invalid.";
    exit(-1);
  }
  try {
    doc["image"] >> mapfname;
    // TODO: make this path-handling more robust
    if (mapfname.size() == 0) {
      std::cout << "The image tag cannot be an empty string."<<std::endl;
      exit(-1);
    }
    /* qiao@2016.06.25
    // [TODO] use relative path...
    if (mapfname[0] != '/') {
      // dirname can modify what you pass it
      char* fname_copy = strdup(fname.c_str());
      mapfname = std::string(dirname(fname_copy)) + '/' + mapfname;
      free(fname_copy);
    }
    */
  } catch (std::exception e) {
       std::cout << "The map does not contain an image tag or it is invalid.";
       exit(-1);
  }

  return LoadMapFromFile(mapfname.c_str(), res, negate, occ_th, free_th, origin, mode, gridmap);
}

bool
LoadMapFromFile(const char* fname, double res, bool negate, double occ_th, double free_th, double* origin, MapMode mode, std::shared_ptr<GridMap> map) {
  cv::Mat img;

  // Load the image using SDL.  If we get NULL back, the image load failed.
  img = cv::imread(fname);
  if(img.empty()) {
    std::string errmsg = std::string("failed to open image file \"") + std::string(fname) + std::string("\"");
    throw std::runtime_error(errmsg);
  }
  // Copy the image data into the map structure
  map.reset();
  map = std::make_shared<GridMap>();
  map->width = img.cols;
  map->height = img.rows;
  map->resolution = res;
  map->origin_x = *(origin);
  map->origin_y = *(origin+1);
  map->origin_th = *(origin+2);

  // Allocate space to hold the data
  map->data.resize(map->width * map->height);

  // Get values that we'll need to iterate through the pixels
  //rowstride = img->pitch;
  int avg_channels;
  // NOTE: Trinary mode still overrides here to preserve existing behavior.
  // Alpha will be averaged in with color channels when using trinary mode.
  if(mode == MapMode::TRINARY || img.channels() == 3 || img.channels() == 1) {
    avg_channels = img.channels();
  } else if(img.channels() == 4) // RGB+alpha
  {
    avg_channels = 3;
  } else {
    std::cout << "LoadMapFromFile() finds an unexpected number of channels";
    return false;
  }

  // Copy pixel data into the map structure
  for(int j = 0; j < (int)map->height; ++j) {
    for(int i = 0; i < (int)map->width; ++i) {
      // Compute mean of RGB for this pixel
      int color_sum = 0;
      int alpha = 1;
      if(img.channels() == 4) {
        color_sum = img.at<cv::Vec4b>(j, i)[0] + img.at<cv::Vec4b>(j, i)[1] + img.at<cv::Vec4b>(j, i)[2];
        alpha = img.at<cv::Vec4b>(j, i)[3];
      } else if(img.channels() == 3) {
        color_sum = img.at<cv::Vec3b>(j, i)[0] + img.at<cv::Vec3b>(j, i)[1] + img.at<cv::Vec3b>(j, i)[2];
      } else {
        color_sum = img.at<unsigned char>(j, i);
      }
      double color_avg = color_sum / (double)avg_channels;

      if(negate) color_avg = 255 - color_avg;

      unsigned char value;
      if(mode == MapMode::RAW) {
        value = color_avg;
        img.at<unsigned char>(j, i) = value;
        continue;
      }

      // If negate is true, we consider blacker pixels free, and whiter
      // pixels free.  Otherwise, it's vice versa.
      const double occ = (255 - color_avg) / 255.0;

      // Apply thresholds to RGB means to determine occupancy values for
      // map.  Note that we invert the graphics-ordering of the pixels to
      // produce a map with cell (0,0) in the lower-left corner.
      if(occ > occ_th) {
        value = 100;
      } else if(occ < free_th) {
        value = 0;
      } else if(mode == MapMode::TRINARY || alpha < 1.0) {
        value = -1;
      } else {
        double ratio = (occ - free_th) / (occ_th - free_th);
        value = 99 * ratio;
      }
      map->data[MAP_IDX(map->width,i,map->height - j - 1)] = value;
      // debug only
      std::cout << (unsigned int)value << " ";
    }
    std::cout << std::endl;
  }
  TRACE_FUNC
  return true;
}

} // namespace map_server
