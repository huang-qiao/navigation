#include "amcl_core.hpp"
#include "util.hpp"

#include <cassert>

#define NEW_UNIFORM_SAMPLING 1

std::vector<std::pair<int, int>> AmclCore::free_space_indices;

Pose AmclCore::UniformPoseGenerator(void *arg) {
  // MapPtr map = (MapPtr)arg;
  Map *map_raw = (Map *)arg;
  MapPtr map = std::make_shared<Map>(*map_raw);
#if NEW_UNIFORM_SAMPLING
  std::cout << "free space indices = " << free_space_indices.size()
            << std::endl;
  unsigned int rand_index = drand48() * free_space_indices.size();
  std::pair<int, int> free_point = free_space_indices[rand_index];
  Pose p;
  p.v[0] = map->toWorldX(free_point.first);
  p.v[1] = map->toWorldY(free_point.second);
  p.v[2] = drand48() * 2 * M_PI - M_PI;
#else
  double min_x, max_x, min_y, max_y;

  min_x = (map->size_x * map->scale) / 2.0 - map->origin_x;
  max_x = (map->size_x * map->scale) / 2.0 + map->origin_x;
  min_y = (map->size_y * map->scale) / 2.0 - map->origin_y;
  max_y = (map->size_y * map->scale) / 2.0 + map->origin_y;

  Pose p;

  ROS_DEBUG("Generating new uniform sample");
  for (;;) {
    p.v[0] = min_x + drand48() * (max_x - min_x);
    p.v[1] = min_y + drand48() * (max_y - min_y);
    p.v[2] = drand48() * 2 * M_PI - M_PI;
    // Check that it's a free cell
    int i, j;
    i = MAP_GXWX(map, p.v[0]);
    j = MAP_GYWY(map, p.v[1]);
    if (MAP_VALID(map, i, j) &&
        (map->cells[map->toIndex(i, j)].occ_state == -1))
      break;
  }
#endif
  return p;
}

AmclCore::AmclCore(double d_thresh, double a_thresh, size_t resample_interval)
    : d_thresh_(d_thresh),
      a_thresh_(a_thresh),
      resample_interval_(resample_interval),
      resample_count_(0) {
  get_first_odom_ = false;
  odom_update_ = false;
  lasers_update_.clear();
  force_update_ = false;
  force_publish_ = false;

  // set default value
  laser_model_ = amcl::LaserModel::LIKELIHOOD_FIELD;
  odom_model_ = amcl::OdomModel::DIFF;

  // debug: display params
  std::cout << "d_thresh_ = " << d_thresh_ << std::endl;
  std::cout << "a_thresh_ = " << a_thresh_ << std::endl;
  std::cout << "resample_interval_ = " << resample_interval_ << std::endl;
}

bool AmclCore::setLaserModel(std::string laser_model_type) {
  if (laser_model_type == "beam") {
    laser_model_ = amcl::LaserModel::BEAM;
  } else if (laser_model_type == "likelihood_field") {
    laser_model_ = amcl::LaserModel::LIKELIHOOD_FIELD;
  } else if (laser_model_type == "likelihood_field_prob") {
    laser_model_ = amcl::LaserModel::LIKELIHOOD_FIELD_PROB;
  } else {
    std::cout << "Unknown laser model type \"" << laser_model_type
              << "\"; defaulting to likelihood_field model" << std::endl;
    laser_model_ = amcl::LaserModel::LIKELIHOOD_FIELD;
  }

  return true;
}

bool AmclCore::setOdomModel(std::string odom_model_type) {
  if (odom_model_type == "diff") {
    odom_model_ = amcl::OdomModel::DIFF;
  } else if (odom_model_type == "omni") {
    odom_model_ = amcl::OdomModel::OMNI;
  } else if (odom_model_type == "diff-corrected") {
    odom_model_ = amcl::OdomModel::DIFF_CORRECTED;
  } else if (odom_model_type == "omni-corrected") {
    odom_model_ = amcl::OdomModel::OMNI_CORRECTED;
  } else {
    std::cout << "Unknown odom model type \"" << odom_model_type
              << "\"; defaulting to diff model" << std::endl;
    odom_model_ = amcl::OdomModel::DIFF;
  }

  return true;
}

bool AmclCore::initParticleFilter(int min_samples, int max_samples,
                                  double alpha_slow, double alpha_fast,
                                  double pop_err, double pop_z,
                                  Pose init_pose_mean,
                                  Covariance init_pose_cov) {
  std::cout << "min_samples: " << min_samples << std::endl;
  std::cout << "max_samples: " << max_samples << std::endl;
  std::cout << "alpha_slow: " << alpha_slow << std::endl;
  std::cout << "alpha_fast: " << alpha_fast << std::endl;
  std::cout << "pop_err: " << pop_err << std::endl;
  std::cout << "pop_z: " << pop_z << std::endl;
  std::cout << "init pose: " << init_pose_mean.v[0] << ","
            << init_pose_mean.v[1] << "," << init_pose_mean.v[2] << std::endl;
#if NEW_UNIFORM_SAMPLING
  // Index of free space
  free_space_indices.resize(0);
  for (int i = 0; i < map_->size_x; i++) {
    for (int j = 0; j < map_->size_y; j++) {
      if (map_->cells[map_->toIndex(i, j)]->occ_state == -1) {
        free_space_indices.push_back(std::make_pair(i, j));
      }
    }
  }
#endif
  pf_ = new ParticleFilter(min_samples, max_samples, alpha_slow, alpha_fast,
                           (pf_init_model_fn_t)AmclCore::UniformPoseGenerator,
                           (void *)map_.get());
  pf_->pop_err_ = pop_err;
  pf_->pop_z_ = pop_z;
  pf_->init(init_pose_mean, init_pose_cov);
  get_first_odom_ = false;
  odom_update_ = false;
  lasers_update_.clear();

  return true;
}

void AmclCore::registerCallbacks(Callbacks *cb) { cbs_.push_back(cb); }

bool AmclCore::initOdom(double alpha1, double alpha2, double alpha3,
                        double alpha4, double alpha5) {
  odom_.reset();
  odom_ = std::make_shared<amcl::AMCLOdom>();
  odom_->SetModel(odom_model_, alpha1, alpha2, alpha3, alpha4, alpha5);

  return true;
}

bool AmclCore::initLaser(size_t max_beams, double z_hit, double z_short,
                         double z_max, double z_rand, double sigma_hit,
                         double lambda_short, double chi_outlier,
                         double max_occ_dist, bool do_beamskip,
                         double beam_skip_distance, double beam_skip_threshold,
                         double beam_skip_error_threshold) {
  if (!map_) {
    std::cerr << "initLaser: map_ not set!";
    return false;
  }

  laser_.reset();  // delete laser_;
  laser_ = std::make_shared<amcl::AMCLLaser>(max_beams, map_);
  if (laser_model_ == amcl::LaserModel::BEAM)
    laser_->SetModelBeam(z_hit, z_short, z_max, z_rand, sigma_hit, lambda_short,
                         chi_outlier);
  else if (laser_model_ == amcl::LaserModel::LIKELIHOOD_FIELD_PROB) {
    std::cout << "Initializing likelihood field model; this can take some time "
                 "on large maps..."
              << std::endl;
    laser_->SetModelLikelihoodFieldProb(
        z_hit, z_rand, sigma_hit, max_occ_dist, do_beamskip, beam_skip_distance,
        beam_skip_threshold, beam_skip_error_threshold);
    std::cout << "Done initializing likelihood field model with probabilities."
              << std::endl;
  } else if (laser_model_ == amcl::LaserModel::LIKELIHOOD_FIELD) {
    std::cout << "Initializing likelihood field model; this can take some time "
                 "on large maps..."
              << std::endl;
    laser_->SetModelLikelihoodField(z_hit, z_rand, sigma_hit, max_occ_dist);
    std::cout << "Done initializing likelihood field model." << std::endl;
  }

  return true;
}

bool AmclCore::initMap(size_t width, size_t height, double scale, double orig_x,
                       double orig_y, std::vector<int> map_data) {
  assert(width * height == map_data.size());

  map_.reset();
  map_ = CreateMap();
  map_->size_x = width;
  map_->size_y = height;
  map_->scale = scale;
  map_->origin_x = orig_x;
  map_->origin_y = orig_y;
  map_->cells.resize(width * height);
  std::fill(map_->cells.begin(), map_->cells.end(), std::make_shared<Cell>());
  for (size_t i = 0; i < map_data.size(); i++) {
    if (map_data[i] == 0) {
      map_->cells[i]->occ_state = -1;
    } else if (map_data[i] == 100) {
      map_->cells[i]->occ_state = +1;
    } else {
      map_->cells[i]->occ_state = 0;
    }
  }

  return true;
}

amcl::AMCLLaserPtr AmclCore::getLaser(std::string laser_frame_id) {
  //
  auto it = frame_to_laser_.find(laser_frame_id);
  if (it == frame_to_laser_.end()) {
    std::cerr << "getLaser: " << laser_frame_id << " not found" << std::endl;
    //
    return nullptr;
  }
  int laser_index = it->second;
  // std::cout << "getLaser: laser_index = " << laser_index << std::endl;
  return lasers_[laser_index];
}

bool AmclCore::registerLaser(std::string laser_frame_id, Pose laser_pose) {
  if (getLaser(laser_frame_id)) {
    // std::cout << laser_frame_id << "is already registered" << std::endl;
    return false;
  }
  lasers_.push_back(std::make_shared<amcl::AMCLLaser>(*laser_));
  lasers_update_.push_back(true);
  size_t laser_index = frame_to_laser_.size();
  // std::cout << "registerLaser: laser_index = " << laser_index << std::endl;
  lasers_[laser_index]->SetLaserPose(laser_pose);
  lasers_[laser_index]->SetLaserFrameId(laser_frame_id);
  frame_to_laser_[laser_frame_id] = laser_index;

  return true;
}

void AmclCore::clearLasers() {
  lasers_.clear();
  frame_to_laser_.clear();
  lasers_update_.clear();
  // if (laser_) {
  //  laser_.reset();
  //}
}
void AmclCore::clearMap() {
  if (map_) {
    map_.reset();
  }
}
void AmclCore::clearOdom() {
  if (odom_) {
    odom_.reset();
  }
}
void AmclCore::clearParticleFilter() {
  if (pf_) {
    delete pf_;
    pf_ = NULL;
    // pf_.reset();
  }
}

bool AmclCore::pfSetUniform() { return true; }

bool AmclCore::pfSetInitPose(Pose pf_pose_mean, Covariance pf_pose_cov) {
  pf_->init(pf_pose_mean, pf_pose_cov);
  get_first_odom_ = false;

  return true;
}

bool AmclCore::pfUpdateOdom(Pose pose) {
  if (!map_) {
    std::cerr << "pfUpdateOdom: map_ not init yet!" << std::endl;
    return false;
  }

  Pose delta = Pose(0.0, 0.0, 0.0);
  if (get_first_odom_) {
    // compute change in pose
    delta.v[0] = pose.v[0] - latest_odom_pose_.v[0];
    delta.v[1] = pose.v[1] - latest_odom_pose_.v[1];
    delta.v[2] = angle_diff(pose.v[2], latest_odom_pose_.v[2]);

    // See if we should update the filter
    bool update = fabs(delta.v[0]) > d_thresh_ ||
                  fabs(delta.v[1]) > d_thresh_ || fabs(delta.v[2]) > a_thresh_;
    update = update || force_update_;
    force_update_ = false;

    // Set the laser update flags
    if (update) {
      // std::cout << "set sensor update to TRUE" << std::endl;
      odom_update_ = true;
      std::fill(lasers_update_.begin(), lasers_update_.end(), true);
    }
  }

  force_publish_ = false;
  if (!get_first_odom_) {
    std::cout << "initializing pf, set pf_odom_pose" << std::endl;
    // Pose at last filter update
    latest_odom_pose_ = pose;

    // Filter is now initialized
    get_first_odom_ = true;

    // force publish core information
    force_publish_ = true;

    resample_count_ = 0;
  } else {
    if (odom_update_) {
      amcl::AMCLOdomDataPtr odata = std::make_shared<amcl::AMCLOdomData>();
      odata->pose = pose;
      odata->delta = delta;
      for (auto cb : cbs_) {
        cb->onUpdateOdom(odata);
      }
      odom_->UpdateAction(
          pf_, std::dynamic_pointer_cast<amcl::AMCLSensorData>(odata));
      odom_update_ = false;
      // update latest odom
      latest_odom_pose_ = pose;
    }
  }
  // [TODO] design a resample algorithm after update odometry...

  if (!force_update_) {
    SampleSet *set = pf_->sets + pf_->current_set_;
    for (auto cb : cbs_) {
      cb->onUpdateParticles(set);
    }
  }

  return true;
}

bool AmclCore::pfUpdateLaser(amcl::AMCLLaserDataPtr ldata) {
  if (!map_) {
    std::cerr << "pfUpdateLaser: map_ not init yet!" << std::endl;

    return false;
  }
  amcl::AMCLLaserPtr laser =
      std::dynamic_pointer_cast<amcl::AMCLLaser>(ldata->sensor);
  std::string laser_frame_id = laser->GetLaserFrameId();
  auto it = frame_to_laser_.find(laser_frame_id);
  if (it == frame_to_laser_.end()) {
    std::cout << "pfUpdateLaser: " << laser_frame_id << " not found"
              << std::endl;
    return false;
  }
  int laser_index = it->second;
  // std::cout << "laser_index = " << laser_index << std::endl;
  // std::cout << "lasers_update_.size = " << lasers_update_.size() << std::endl;
  if (!lasers_update_[laser_index]) {
    // std::cout << "pfUpdateLaser: laser update flag not set" << std::endl;
    return false;
  }

  laser->UpdateSensor(pf_,
                      std::dynamic_pointer_cast<amcl::AMCLSensorData>(ldata));
  lasers_update_[laser_index] = false;

  bool resampled = false;
  if (!(++resample_count_ % resample_interval_)) {
    pf_->updateResample();
    resampled = true;
  }

  if (resampled || force_publish_) {
    SampleSet *set = pf_->sets + pf_->current_set_;
    for (auto cb : cbs_) {
      cb->onUpdateParticles(set);
    }
  }

  return true;
}
