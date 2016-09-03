#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "map.hpp"


////////////////////////////////////////////////////////////////////////////
// Load an occupancy grid
int Map::loadOcc(const char *filename, double scale, int negate)
{
  FILE *file;
  char magic[3];
  int i, j;
  int ch, occ;
  int width, height, depth;
  CellPtr cell;

  // Open file
  file = fopen(filename, "r");
  if (file == NULL)
  {
    fprintf(stderr, "%s: %s\n", strerror(errno), filename);
    return -1;
  }

  // Read ppm header

  if ((fscanf(file, "%10s \n", magic) != 1) || (strcmp(magic, "P5") != 0))
  {
    fprintf(stderr, "incorrect image format; must be PGM/binary");
    return -1;
  }

  // Ignore comments
  while ((ch = fgetc(file)) == '#')
    while (fgetc(file) != '\n');
  ungetc(ch, file);

  // Read image dimensions
  if(fscanf(file, " %d %d \n %d \n", &width, &height, &depth) != 3)
  {
    fprintf(stderr, "Failed ot read image dimensions");
    return -1;
  }

  // Allocate space in the map
  if (cells.empty())
  {
    this->scale = scale;
    this->size_x = width;
    this->size_y = height;
    //map->cells = (MapCell*)calloc(width * height, sizeof(map->cells[0]));
    this->cells.resize(width * height);
    std::fill(this->cells.begin(), this->cells.end(), std::make_shared<Cell>());
  }
  else
  {
    if (width != this->size_x || height != this->size_y)
    {
      //PLAYER_ERROR("map dimensions are inconsistent with prior map dimensions");
      return -1;
    }
  }

  // Read in the image
  for (j = height - 1; j >= 0; j--)
  {
    for (i = 0; i < width; i++)
    {
      ch = fgetc(file);

      // Black-on-white images
      if (!negate)
      {
        if (ch < depth / 4)
          occ = +1;
        else if (ch > 3 * depth / 4)
          occ = -1;
        else
          occ = 0;
      }

      // White-on-black images
      else
      {
        if (ch < depth / 4)
          occ = -1;
        else if (ch > 3 * depth / 4)
          occ = +1;
        else
          occ = 0;
      }

      if (!isValid(i, j))
        continue;
      //cell = map->cells + MAP_INDEX(map, i, j);
      cell = cells[toIndex(i, j)];
      cell->occ_state = occ;
    }
  }

  fclose(file);

  return 0;
}


////////////////////////////////////////////////////////////////////////////
// Load a wifi signal strength map
/*
int map_load_wifi(MapPtr map, const char *filename, int index)
{
  FILE *file;
  char magic[3];
  int i, j;
  int ch, level;
  int width, height, depth;
  CellPtr cell;

  // Open file
  file = fopen(filename, "r");
  if (file == NULL)
  {
    fprintf(stderr, "%s: %s\n", strerror(errno), filename);
    return -1;
  }

  // Read ppm header
  fscanf(file, "%10s \n", magic);
  if (strcmp(magic, "P5") != 0)
  {
    fprintf(stderr, "incorrect image format; must be PGM/binary");
    return -1;
  }

  // Ignore comments
  while ((ch = fgetc(file)) == '#')
    while (fgetc(file) != '\n');
  ungetc(ch, file);

  // Read image dimensions
  fscanf(file, " %d %d \n %d \n", &width, &height, &depth);

  // Allocate space in the map
  if (map->cells == NULL)
  {
    map->size_x = width;
    map->size_y = height;
    map->cells = calloc(width * height, sizeof(map->cells[0]));
  }
  else
  {
    if (width != map->size_x || height != map->size_y)
    {
      //PLAYER_ERROR("map dimensions are inconsistent with prior map dimensions");
      return -1;
    }
  }

  // Read in the image
  for (j = height - 1; j >= 0; j--)
  {
    for (i = 0; i < width; i++)
    {
      ch = fgetc(file);

      if (!MAP_VALID(map, i, j))
        continue;

      if (ch == 0)
        level = 0;
      else
        level = ch * 100 / 255 - 100;

      cell = map->cells + MAP_INDEX(map, i, j);
      cell->wifi_levels[index] = level;
    }
  }

  fclose(file);

  return 0;
}
*/



