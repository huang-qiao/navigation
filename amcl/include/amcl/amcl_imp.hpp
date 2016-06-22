#pragma once

#include "amcl/pf/pf.h"
#include "amcl/map/map.h"
#include "amcl/sensors/amcl_laser.h"
#include "sensors/amcl_sensor.h"
#include "sensors/amcl_odom.h"

#include "Eigen/Geometry"

#include <cmath>
#include <map>
#include <vector>
#include <iostream>

#define NEW_UNIFORM_SAMPLING 1

namespace amcl {

// Pose hypothesis
typedef struct
{
  // Total weight (weights sum to 1)
  double weight;

  // Mean of pose esimate
  pf_vector_t pf_pose_mean;

  // Covariance of pose estimate
  pf_matrix_t pf_pose_cov;

} amcl_hyp_t;

double
normalize(double z)
{
  return atan2(sin(z),cos(z));
}

double
angle_diff(double a, double b)
{
  double d1, d2;
  a = normalize(a);
  b = normalize(b);
  d1 = a-b;
  d2 = 2*M_PI - fabs(d1);
  if(d1 > 0)
    d2 *= -1.0;
  if(fabs(d1) < fabs(d2))
    return(d1);
  else
    return(d2);
}

class AmclImp;
class AmclImpCallback {
public:
  virtual void onUpdatePointCloud(pf_sample_set_t* set) = 0;
  virtual void onUpdatePose(pf_vector_t pf_pose_mean, pf_matrix_t cov, double time_sec) = 0;
  virtual void onUpdateTf(pf_vector_t pf_pose_mean, double time_sec) = 0;
  virtual void onQueryInitPose(AmclImp* amcl) = 0;
  virtual void onSaveCurrentPose(pf_vector_t pf_pose, pf_matrix_t pf_cov) = 0;
};

class AmclImp
{
public:
  AmclImp();
  ~AmclImp(){}

  inline map_t* getMap() const {
    std::cout << __PRETTY_FUNCTION__ << "++" << std::endl;
    return map_;
  }

  inline void setInitPose(double* init_pose, double* init_cov)
  {
    std::cout << __PRETTY_FUNCTION__ << "++" << std::endl;
    init_pose_[0] = init_pose[0];
    init_pose_[1] = init_pose[1];
    init_pose_[2] = init_pose[2];
    init_cov_[0] = init_cov[0];
    init_cov_[1] = init_cov[1];
    init_cov_[2] = init_cov[2];
  }

  inline void initPoseHypothesis(pf_vector_t pf_init_pose_mean, pf_matrix_t pf_init_pose_cov)
  {
    std::cout << __PRETTY_FUNCTION__ << "++" << std::endl;
    delete initial_pose_hyp_;
    initial_pose_hyp_ = new amcl_hyp_t();
    initial_pose_hyp_->pf_pose_mean = pf_init_pose_mean;
    initial_pose_hyp_->pf_pose_cov = pf_init_pose_cov;
    applyInitialPose();
  }

  inline void freeMapDependentMemory()
  {
    std::cout << __PRETTY_FUNCTION__ << "++" << std::endl;
    if( map_ != NULL ) {
      map_free( map_ );
      map_ = NULL;
    }
    if( pf_ != NULL ) {
      pf_free( pf_ );
      pf_ = NULL;
    }

    delete odom_;
    odom_ = NULL;
    delete laser_;
    laser_ = NULL;
  }

  inline void initMap(map_t* map)
  {
    lasers_.clear();
    lasers_update_.clear();
    frame_to_laser_.clear();

    map_ = map;

  #if NEW_UNIFORM_SAMPLING
    // Index of free space
    free_space_indices.resize(0);
    for(int i = 0; i < map_->size_x; i++)
      for(int j = 0; j < map_->size_y; j++)
        if(map_->cells[MAP_INDEX(map_,i,j)].occ_state == -1)
          free_space_indices.push_back(std::make_pair(i,j));
  #endif
    // Create the particle filter
    pf_ = pf_alloc(min_particles_, max_particles_,
                   alpha_slow_, alpha_fast_,
                   (pf_init_model_fn_t)AmclImp::uniformPoseGenerator,
                   (void *)map_);
    pf_->pop_err = pf_err_;
    pf_->pop_z = pf_z_;

    // Initialize the filter
    //updatePoseFromServer();
    queryInitPose();
    pf_vector_t pf_init_pose_mean = pf_vector_zero();
    pf_init_pose_mean.v[0] = init_pose_[0];
    pf_init_pose_mean.v[1] = init_pose_[1];
    pf_init_pose_mean.v[2] = init_pose_[2];
    pf_matrix_t pf_init_pose_cov = pf_matrix_zero();
    pf_init_pose_cov.m[0][0] = init_cov_[0];
    pf_init_pose_cov.m[1][1] = init_cov_[1];
    pf_init_pose_cov.m[2][2] = init_cov_[2];
    pf_init(pf_, pf_init_pose_mean, pf_init_pose_cov);
    pf_init_ = false;

    // Instantiate the sensor objects
    // Odometry
    delete odom_;
    odom_ = new AMCLOdom();
    //ROS_ASSERT(odom_);
    odom_->SetModel( odom_model_type_, alpha1_, alpha2_, alpha3_, alpha4_, alpha5_ );
    // Laser
    delete laser_;
    laser_ = new AMCLLaser(max_beams_, map_);
    //ROS_ASSERT(laser_);
    if(laser_model_type_ == LASER_MODEL_BEAM)
      laser_->SetModelBeam(z_hit_, z_short_, z_max_, z_rand_,
                           sigma_hit_, lambda_short_, 0.0);
    else if(laser_model_type_ == LASER_MODEL_LIKELIHOOD_FIELD_PROB){
      printf("Initializing likelihood field model; this can take some time on large maps...");
      laser_->SetModelLikelihoodFieldProb(z_hit_, z_rand_, sigma_hit_,
                                          laser_likelihood_max_dist_,
                                          do_beamskip_, beam_skip_distance_,
                                          beam_skip_threshold_, beam_skip_error_threshold_);
      printf("Done initializing likelihood field model.");
    }
    else
    {
      printf("Initializing likelihood field model; this can take some time on large maps...");
      laser_->SetModelLikelihoodField(z_hit_, z_rand_, sigma_hit_,
                                      laser_likelihood_max_dist_);
      printf("Done initializing likelihood field model.");
    }

    // In case the initial pose message arrived before the first map,
    // try to apply the initial pose now that the map has arrived.
    applyInitialPose();
    first_map_received_ = true;
  }

  /**
   * If initial_pose_hyp_ and map_ are both non-null, apply the initial
   * pose to the particle filter state.  initial_pose_hyp_ is deleted
   * and set to NULL after it is used.
   */
  inline void applyInitialPose()
  {
    std::cout << __PRETTY_FUNCTION__ << "++" << std::endl;
    if( initial_pose_hyp_ != NULL && map_ != NULL ) {
      pf_init(pf_, initial_pose_hyp_->pf_pose_mean, initial_pose_hyp_->pf_pose_cov);
      pf_init_ = false;

      delete initial_pose_hyp_;
      initial_pose_hyp_ = NULL;
    }
  }

  inline bool isLaserUpdate(std::string frame_id)
  {
    std::cout << __PRETTY_FUNCTION__ << "++" << std::endl;
    int laser_index = frame_to_laser_[frame_id];
    return lasers_update_[laser_index];
  }

  inline pf_vector_t getLatestOdomPose() const {
    std::cout << __PRETTY_FUNCTION__ << "++" << std::endl;
    return latest_pose_;
  }

  inline AMCLLaser* getLaser(std::string frame_id)
  {
    std::cout << __PRETTY_FUNCTION__ << "++" << std::endl;
    int laser_index = frame_to_laser_[frame_id];
    return lasers_[laser_index];
  }

  inline double getLaserMaxRange() {
    std::cout << __PRETTY_FUNCTION__ << "++" << std::endl;
    return laser_max_range_;
  }
  inline double getLaserMinRange() {
    std::cout << __PRETTY_FUNCTION__ << "++" << std::endl;
    return laser_min_range_;
  }

  inline bool doGlobalLocalization()
  {
    std::cout << __PRETTY_FUNCTION__ << "++" << std::endl;
    printf("Initializing with uniform distribution");
    pf_init_model(pf_, (pf_init_model_fn_t)AmclImp::uniformPoseGenerator, (void *)map_);
    printf("Global initialisation done!");
    pf_init_ = false;

    return true;
  }

  inline bool sentFirstTransform() {
    std::cout << __PRETTY_FUNCTION__ << "++" << std::endl;
    return sent_first_transform_;
  }

  inline bool doNoMotionUpdate()
  {
    std::cout << __PRETTY_FUNCTION__ << "++" << std::endl;
    m_force_update = true;
    printf("Requesting no-motion update");
    return true;
  }

  // callback
  inline void onUpdatePointCloud(pf_sample_set_t* set) {
    std::cout << __PRETTY_FUNCTION__ << "++" << std::endl;
    for (auto cb : cb_) {
      cb->onUpdatePointCloud(set);
    }
  }

  inline void onUpdatePose(pf_vector_t pf_pose_mean, pf_matrix_t cov, double time_sec)
  {
    std::cout << __PRETTY_FUNCTION__ << "++" << std::endl;
    for (auto cb : cb_) {
      cb->onUpdatePose(pf_pose_mean, cov, time_sec);
    }
  }

  inline void onUpdateTf(pf_vector_t pf_pose_mean, double time_sec)
  {
    std::cout << __PRETTY_FUNCTION__ << "++" << std::endl;
    for (auto cb : cb_) {
      cb->onUpdateTf(pf_pose_mean, time_sec);
    }
    sent_first_transform_ = true;
  }

  inline void queryInitPose()
  {
    std::cout << __PRETTY_FUNCTION__ << "++" << std::endl;
    for (auto cb : cb_) {
      cb->onQueryInitPose(this);
    }
  }

  inline void onSavePose()
  {
    std::cout << __PRETTY_FUNCTION__ << "++" << std::endl;
    for (auto cb : cb_) {
      cb->onSaveCurrentPose(latest_pose_, latest_cov_);
    }
  }

  inline void registerCallback(AmclImpCallback* cb)
  {
    std::cout << __PRETTY_FUNCTION__ << "++" << std::endl;
    cb_.push_back(cb);
  }

  // predict
  void processSensorPose(std::string frame_id, Eigen::Isometry3d laser_pose, pf_vector_t odom_pose);

  // adjust
  void processLaserData(std::string frame_id, AMCLLaserData* ldata, pf_vector_t odom_pose);
  // Pose-generating function used to uniformly distribute particles over
  // the map
  static pf_vector_t uniformPoseGenerator(void* arg);
#if NEW_UNIFORM_SAMPLING
  static std::vector<std::pair<int,int> > free_space_indices;
#endif

private:
  bool sent_first_transform_;

  map_t* map_;
  char* mapdata;
  int sx, sy;
  double resolution;

  std::vector< AMCLLaser* > lasers_;
  std::vector< bool > lasers_update_;
  std::map< std::string, int > frame_to_laser_;

  // Particle filter
  pf_t *pf_;
  double pf_err_, pf_z_;
  bool pf_init_;
  pf_vector_t pf_odom_pose_;
  double d_thresh_, a_thresh_;
  int resample_interval_;
  int resample_count_;
  double laser_min_range_;
  double laser_max_range_;

  // Nomotion update control
  // used to temporarily let amcl update samples even when no motion occurs...
  bool m_force_update;

  AMCLOdom* odom_;
  AMCLLaser* laser_;

  // callbacks
  std::vector<AmclImpCallback*> cb_;

  amcl_hyp_t* initial_pose_hyp_;
  bool first_map_only_;
  bool first_map_received_;
  bool first_reconfigure_call_;

  bool latest_pose_valid_;
  pf_vector_t latest_pose_;
  pf_matrix_t latest_cov_;

  int max_beams_, min_particles_, max_particles_;
  double alpha1_, alpha2_, alpha3_, alpha4_, alpha5_;
  double alpha_slow_, alpha_fast_;
  double z_hit_, z_short_, z_max_, z_rand_, sigma_hit_, lambda_short_;
//beam skip related params
  bool do_beamskip_;
  double beam_skip_distance_, beam_skip_threshold_, beam_skip_error_threshold_;
  double laser_likelihood_max_dist_;
  odom_model_t odom_model_type_;
  double init_pose_[3];
  double init_cov_[3];
  laser_model_t laser_model_type_;
  bool force_publication_;
}; // AmclImp

} // namespace amcl
