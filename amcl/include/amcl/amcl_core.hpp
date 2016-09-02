#pragma once

#include "map/map.hpp"
#include "pf/pf.hpp"
#include "sensors/amcl_laser.hpp"
#include "sensors/amcl_odom.hpp"

#include <map>
#include <string>
#include <vector>

#define NEW_UNIFORM_SAMPLING 1

class AmclCore {
 public:
  class Callbacks {
   public:
    virtual void onUpdateOdom(amcl::AMCLOdomDataPtr odata) = 0;
    virtual void onUpdateLaser(amcl::AMCLLaserDataPtr ldata) = 0;
    virtual void onUpdateParticles(SampleSet* set) = 0;
  };

 public:
  AmclCore (double d_thresh, double a_thresh, size_t resample_interval);
  bool setLaserModel(std::string laser_model_type);
  inline void setLaserModel(amcl::LaserModel laser_model) { laser_model_ = laser_model; }
  bool initLaser(size_t max_beams, double z_hit, double z_short, double z_max,
                 double z_rand, double sigma_hit, double lambda_short,
                 double chi_outlier, double max_occ_dist, bool do_beamskip,
                 double beam_skip_distance, double beam_skip_threshold,
                 double beam_skip_error_threshold);
  amcl::AMCLLaserPtr getLaser(std::string laser_frame_id);
  bool registerLaser(std::string laser_frame_id, Pose laser_pose);
  bool setOdomModel(std::string odom_model_type);
  inline void setOdomModel(amcl::OdomModel odom_model) { odom_model_ = odom_model; }
  bool initOdom(double a1, double a2, double a3, double a4, double a5);
  bool initMap(size_t width, size_t height, double scale, double orig_x,
               double orig_y, std::vector<int> map_data);
  bool initParticleFilter(int min_samples_, int max_samples_,
                          double alpha_slow_, double alpha_fast_,
                          double pop_err, double pop_z, Pose init_pose_mean,
                          Covariance init_pose_cov);
  bool pfSetUniform();
  bool pfSetInitPose(Pose pf_pose_mean, Covariance pf_pose_cov);
  inline MapPtr getMap() { return map_; }

  bool pfUpdateOdom(Pose pose);
  bool pfUpdateLaser(amcl::AMCLLaserDataPtr ldata);

  void registerCallbacks(Callbacks* cb);
  inline void clearCallbacks() { cbs_.clear(); }

  void clearMap();
  void clearLasers();
  void clearOdom();
  void clearParticleFilter();

 public:
  // Pose-generating function used to uniformly distribute particles over
  // the map
  static Pose UniformPoseGenerator(void* map);
#if NEW_UNIFORM_SAMPLING
  static std::vector<std::pair<int, int>> free_space_indices;
#endif

 private:
  std::vector<Callbacks*> cbs_;
  ParticleFilter* pf_;
  bool get_first_odom_;
  Pose latest_odom_pose_;

  amcl::AMCLOdomPtr odom_;
  amcl::AMCLLaserPtr laser_;
  MapPtr map_;

  amcl::OdomModel odom_model_;
  amcl::LaserModel laser_model_;

  std::vector<amcl::AMCLLaserPtr> lasers_;
  std::vector<bool> lasers_update_;
  std::map<std::string, int> frame_to_laser_;

  double d_thresh_, a_thresh_;
  bool odom_update_, force_update_, force_publish_;
  size_t resample_count_, resample_interval_;
};
