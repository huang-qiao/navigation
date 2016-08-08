#include <queue>
#include <cmath>
#include <cstdlib>
#include <string>
#include "map.hpp"

#include <cstring>
#include <iostream>
#define __FILENAME__ \
  (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)
#define TRACE_FUNC                                                     \
  do {                                                                 \
    std::cout << __FILENAME__ << ":" << __LINE__ << " in " << __func__ \
              << std::endl;                                            \
  } while (0);

#define TRACE_FUNC_ENTER                                               \
  do {                                                                 \
    std::cout << __FILENAME__ << ":" << __LINE__ << " in " << __func__ \
              << ": [ENTER]" << std::endl;                             \
  } while (0);

#define TRACE_FUNC_EXIT                                                \
  do {                                                                 \
    std::cout << __FILENAME__ << ":" << __LINE__ << " in " << __func__ \
              << ": [EXIT]" << std::endl;                              \
  } while (0);

struct CellData
{
  amcl::map::MapPtr map_;
  unsigned int i_, j_;
  unsigned int src_i_, src_j_;
};

bool operator<(const CellData& a, const CellData& b)
{
  return a.map_->cells[a.map_->toIndex(a.i_, a.j_)]->occ_dist > a.map_->cells[b.map_->toIndex(b.i_, b.j_)]->occ_dist;
}

class CachedDistanceMap
{
public:
  inline CachedDistanceMap(const double &scale, const double &max_dist) : scale_(scale), max_dist_(max_dist)
  {
    cell_radius_ = max_dist / scale;
    distances_.resize(cell_radius_ + 2);
    for(int i=0; i<=cell_radius_+1; i++) {
      distances_[i].resize(cell_radius_+2);
      for(int j=0; j<=cell_radius_+1; j++) {
        distances_[i][j] = sqrt(i*i + j*j);
      }
    }
  }
  inline ~CachedDistanceMap()
  {
    if(!distances_.empty()){
      for(int i=0; i<=cell_radius_+1; i++) {
        distances_[i].clear();
      }
      distances_.clear();
    }
  }
public: // [FIXME] should be private
  std::vector<std::vector<double>> distances_;
  double scale_;
  double max_dist_;
  int cell_radius_;
};

//using CachedDistanceMapPtr = std::shared_ptr<CachedDistanceMap>;

CachedDistanceMapPtr
GetDistanceMap(double scale, double max_dist)
{
  static CachedDistanceMapPtr cdm = nullptr;

  if(!cdm || (cdm->scale_ != scale) || (cdm->max_dist_ != max_dist))
  {
    if (cdm) {
      cdm.reset();
    }
    cdm = std::make_shared<CachedDistanceMap>(scale, max_dist);
  }

  return cdm;
}


//void enqueue(map_t* map, unsigned int i, unsigned int j, unsigned int src_i, unsigned int src_j,
//	       std::priority_queue<CellData>& Q, CachedDistanceMap* cdm, unsigned char* marked)
void amcl::map::Map::enqueue(const int &i, const int &j, const int &src_i, const int &src_j, std::priority_queue<CellData> &Q, CachedDistanceMapPtr cdm, std::vector<bool> &marked)
{
  //TRACE_FUNC_ENTER
  if(marked[toIndex(i, j)])
    return;

  unsigned int di = std::abs(i - src_i);
  unsigned int dj = std::abs(j - src_j);
  double distance = cdm->distances_[di][dj];

  if(distance > cdm->cell_radius_)
    return;

  cells[toIndex(i,j)]->occ_dist = distance * scale;

  CellData cell;
  cell.map_ = std::shared_ptr<Map>(this);
  cell.i_ = i;
  cell.j_ = j;
  cell.src_i_ = src_i;
  cell.src_j_ = src_j;

  Q.push(cell);

  marked[toIndex(i, j)] = true;
  //TRACE_FUNC_EXIT
}

// Update the cspace distance values
//void map_update_cspace(map_t *map, double max_occ_dist)
void amcl::map::Map::updateCSpace(const double &max_occ_dist)
{
  TRACE_FUNC_ENTER
  std::vector<bool> marked;
  std::priority_queue<CellData> Q;

  marked.resize(size_x * size_y);
  std::fill(marked.begin(), marked.end(), false);

  max_occ_dist_ = max_occ_dist;

  std::cout << "scale = " << scale << ", max_occ_dist = " << max_occ_dist_ << std::endl;
  CachedDistanceMapPtr cdm = GetDistanceMap(scale, max_occ_dist_);

  // Enqueue all the obstacle cells
  CellData cell;
  cell.map_ = std::shared_ptr<Map>(this);
  for(int i=0; i<size_x; i++) {
    cell.src_i_ = cell.i_ = i;
    for(int j=0; j<size_y; j++) {
      if(cells[toIndex(i, j)]->occ_state == +1) {
        cells[toIndex(i, j)]->occ_dist = 0.0;
	cell.src_j_ = cell.j_ = j;
	marked[toIndex(i, j)] = true;
	Q.push(cell);
      } else {
        cells[toIndex(i, j)]->occ_dist = max_occ_dist;
      }
    }
  }

  while(!Q.empty()) {
    CellData current_cell = Q.top();
    if(current_cell.i_ > 0)
      enqueue(current_cell.i_-1, current_cell.j_,
	      current_cell.src_i_, current_cell.src_j_,
	      Q, cdm, marked);
    if(current_cell.j_ > 0)
      enqueue(current_cell.i_, current_cell.j_-1,
	      current_cell.src_i_, current_cell.src_j_,
	      Q, cdm, marked);
    if((int)current_cell.i_ < size_x - 1)
      enqueue(current_cell.i_+1, current_cell.j_,
	      current_cell.src_i_, current_cell.src_j_,
	      Q, cdm, marked);
    if((int)current_cell.j_ < size_y - 1)
      enqueue(current_cell.i_, current_cell.j_+1,
	      current_cell.src_i_, current_cell.src_j_,
	      Q, cdm, marked);
    Q.pop();
    // debug only
    std::cout << "Q.size = " << Q.size() << std::endl;
  }
  TRACE_FUNC_EXIT
}
