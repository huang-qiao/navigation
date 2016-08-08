#include <sys/types.h>  // required by Darwin
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <unistd.h>

#include "amcl_laser.hpp"

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

// floating point comparasion
static bool Equal(const double &a, const double &b)
{
  return std::fabs(a -b) < 1e-6;
}

static bool Equal(const float &a, const float &b)
{
  return Equal((double)a, (double)b);
}

using namespace amcl;

////////////////////////////////////////////////////////////////////////////////
// Default constructor
AMCLLaser::AMCLLaser(const size_t &max_beams, const map::MapPtr map)
    : AMCLSensor(), max_samples_(0), max_obs_(0), temp_obs(NULL) {
  time_ = 0.0;

  max_beams_ = max_beams;
  map_ = map;

  return;
}

AMCLLaser::~AMCLLaser() {
  if (temp_obs) {
    for (int k = 0; k < max_samples_; k++) {
      delete[] temp_obs[k];
    }
    delete[] temp_obs;
  }
}

void AMCLLaser::setModelBeam(double z_hit, double z_short, double z_max,
                             double z_rand, double sigma_hit,
                             double lambda_short, double chi_outlier) {
  TRACE_FUNC_ENTER
  model_type_ = LASER_MODEL_BEAM;
  z_hit_ = z_hit;
  z_short_ = z_short;
  z_max_ = z_max;
  z_rand_ = z_rand;
  sigma_hit_ = sigma_hit;
  lambda_short_ = lambda_short;
  chi_outlier_ = chi_outlier;
  TRACE_FUNC_EXIT
}

void AMCLLaser::setModelLikelihoodField(double z_hit, double z_rand,
                                        double sigma_hit, double max_occ_dist) {
  TRACE_FUNC_ENTER
  model_type_ = LASER_MODEL_LIKELIHOOD_FIELD;
  z_hit_ = z_hit;
  z_rand_ = z_rand;
  sigma_hit_ = sigma_hit;

  map_->updateCSpace(max_occ_dist);
  TRACE_FUNC_EXIT
}

void AMCLLaser::setModelLikelihoodFieldProb(
    double z_hit, double z_rand, double sigma_hit, double max_occ_dist,
    bool do_beamskip, double beam_skip_distance, double beam_skip_threshold,
    double beam_skip_error_threshold) {
  TRACE_FUNC_ENTER
  model_type_ = LASER_MODEL_LIKELIHOOD_FIELD_PROB;
  z_hit_ = z_hit;
  z_rand_ = z_rand;
  sigma_hit_ = sigma_hit;
  do_beamskip_ = do_beamskip;
  beam_skip_distance_ = beam_skip_distance;
  beam_skip_threshold_ = beam_skip_threshold;
  beam_skip_error_threshold_ = beam_skip_error_threshold;
  map_->updateCSpace(max_occ_dist);
  TRACE_FUNC_EXIT
}

////////////////////////////////////////////////////////////////////////////////
// Apply the laser sensor model
bool AMCLLaser::updateSensor(pf::ParticleFilterPtr pf, AMCLLaserDataPtr data) {
  if (max_beams_ < 2) return false;

  double total_weight = 0.0;
  // Apply the laser sensor model
  if (model_type_ == LASER_MODEL_BEAM) {
    //pf->updateSensor(func1, data);
    total_weight = BeamModel(data, pf->getCurrentSet());
  } else if (model_type_ == LASER_MODEL_LIKELIHOOD_FIELD) {
    //pf_update_sensor(pf, (pf_sensor_model_fn_t)LikelihoodFieldModel, data);
    total_weight = LikelihoodFieldModel(data, pf->getCurrentSet());
  } else if (model_type_ == LASER_MODEL_LIKELIHOOD_FIELD_PROB) {
    //pf_update_sensor(pf, (pf_sensor_model_fn_t)LikelihoodFieldModelProb, data);
    total_weight = LikelihoodFieldModelProb(data, pf->getCurrentSet());
  } else {
    //pf_update_sensor(pf, (pf_sensor_model_fn_t)BeamModel, data);
    total_weight = BeamModel(data, pf->getCurrentSet());
  }

  pf->updateSampleWeights(total_weight);

  return true;
}

////////////////////////////////////////////////////////////////////////////////
// Determine the probability for the given pose
double AMCLLaser::BeamModel(AMCLLaserDataPtr data, pf::SampleSetPtr set) {
  //AMCLLaser *self;
  //int i, j, step;
  //double z, pz;
  //double p;
  //double map_range;
  //double obs_range, obs_bearing;
  //double total_weight;
  //pf_sample_t *sample;
  //pf_vector_t pose;

  //self = (AMCLLaser *)data->sensor;
  AMCLLaserPtr self = std::dynamic_pointer_cast<AMCLLaser>(data->sensor);

  double total_weight = 0.0;

  // Compute the sample weights
  for (size_t j = 0; j < set->samples.size(); j++) {
    pf::SamplePtr sample = set->samples[j];
    pf::Pose pose = sample->pose;

    // Take account of the laser pose relative to the robot
    //pose = pf_vector_coord_add(self->laser_pose, pose);
    pose = pf::Pose::CoordAdd(self->laser_pose_, pose);

    double p = 1.0;

    int step = (data->ranges.size() - 1) / (self->max_beams_ - 1);
    for (size_t i = 0; i < data->ranges.size(); i += step) {
      //obs_range = data->ranges[i][0];
      //obs_bearing = data->ranges[i][1];
      const double obs_range = data->ranges[i].range;
      const double obs_bearing = data->ranges[i].angle;

      // Compute the range according to the map
      //map_range = map_calc_range(self->map, pose.v[0], pose.v[1], pose.v[2] + obs_bearing, data->range_max);
      const double map_range = self->map_->calcRange(pose.x(), pose.y(), pose.a(), data->range_max);

      double pz = 0.0;

      // Part 1: good, but noisy, hit
      const double z = obs_range - map_range;
      pz += self->z_hit_ * exp(-(z * z) / (2 * self->sigma_hit_ * self->sigma_hit_));

      // Part 2: short reading from unexpected obstacle (e.g., a person)
      if (z < 0) {
        pz += self->z_short_ * self->lambda_short_ * exp(-self->lambda_short_ * obs_range);
      }

      // Part 3: Failure to detect obstacle, reported as max-range
      if (Equal(obs_range, data->range_max)) {
          pz += self->z_max_ * 1.0;
      }

      // Part 4: Random measurements
      if (obs_range < data->range_max) {
        pz += self->z_rand_ * 1.0 / data->range_max;
      }
      // TODO: outlier rejection for short readings

      assert(pz <= 1.0);
      assert(pz >= 0.0);
      //      p *= pz;
      // here we have an ad-hoc weighting scheme for combining beam probs
      // works well, though...
      p += pz * pz * pz;
    }

    sample->weight *= p;
    total_weight += sample->weight;
  }

  return (total_weight);
}

//double AMCLLaser::LikelihoodFieldModel(AMCLLaserData *data, pf_sample_set_t *set)
double AMCLLaser::LikelihoodFieldModel(AMCLLaserDataPtr data, pf::SampleSetPtr set)
{
  AMCLLaserPtr self = std::dynamic_pointer_cast<AMCLLaser>(data->sensor);

  double total_weight = 0.0;

  // Compute the sample weights
  for (size_t j = 0; j < set->samples.size(); j++) {
    pf::SamplePtr sample = set->samples[j];
    pf::Pose pose = sample->pose;

    // Take account of the laser pose relative to the robot
    pose = pf::Pose::CoordAdd(self->laser_pose_, pose);

    double p = 1.0;

    // Pre-compute a couple of things
    double z_hit_denom = 2 * self->sigma_hit_ * self->sigma_hit_;
    double z_rand_mult = 1.0 / data->range_max;

    int step = (data->ranges.size() - 1) / (self->max_beams_ - 1);

    // Step size must be at least 1
    if (step < 1) step = 1;

    for (size_t i = 0; i < data->ranges.size(); i += step) {
      double obs_range = data->ranges[i].range;
      double obs_bearing = data->ranges[i].angle;

      // This model ignores max range readings
      if (obs_range >= data->range_max) {
        continue;
      }

      // Check for NaN
      if (std::isnan(obs_range)) {
        continue;
      }

      double pz = 0.0;

      pf::Pose hit;
      // Compute the endpoint of the beam
      hit.v[0] = pose.v[0] + obs_range * cos(pose.v[2] + obs_bearing);
      hit.v[1] = pose.v[1] + obs_range * sin(pose.v[2] + obs_bearing);

      // Convert to map grid coords.
      int mi = self->map_->toGridX(hit.x());
      int mj = self->map_->toGridY(hit.y());

      double z = 0.0;
      // Part 1: Get distance from the hit to closest obstacle.
      // Off-map penalized as max distance
      if (!self->map_->isValid(mi, mj)) {
        z = self->map_->max_occ_dist_;
      } else {
        z = self->map_->cells[self->map_->toIndex(mi, mj)]->occ_dist;
      }
      // Gaussian model
      // NOTE: this should have a normalization of 1/(sqrt(2pi)*sigma)
      pz += self->z_hit_ * exp(-(z * z) / z_hit_denom);

      // Part 2: random measurements
      pz += self->z_rand_ * z_rand_mult;

      // TODO: outlier rejection for short readings

      assert(pz <= 1.0);
      assert(pz >= 0.0);
      //      p *= pz;
      // here we have an ad-hoc weighting scheme for combining beam probs
      // works well, though...
      p += pz * pz * pz;
    }

    sample->weight *= p;
    total_weight += sample->weight;
  }

  return (total_weight);
}

//double AMCLLaser::LikelihoodFieldModelProb(AMCLLaserData *data, pf_sample_set_t *set)
double AMCLLaser::LikelihoodFieldModelProb(AMCLLaserDataPtr data, pf::SampleSetPtr set)
{
  AMCLLaserPtr self = std::dynamic_pointer_cast<AMCLLaser>(data->sensor);

  double total_weight = 0.0;

  int step = ceil((data->ranges.size()) / static_cast<double>(self->max_beams_));

  // Step size must be at least 1
  if (step < 1) step = 1;

  // Pre-compute a couple of things
  double z_hit_denom = 2 * self->sigma_hit_ * self->sigma_hit_;
  double z_rand_mult = 1.0 / data->range_max;

  double max_dist_prob =
      exp(-(self->map_->max_occ_dist_ * self->map_->max_occ_dist_) / z_hit_denom);

  // Beam skipping - ignores beams for which a majoirty of particles do not
  // agree with the map
  // prevents correct particles from getting down weighted because of unexpected
  // obstacles
  // such as humans
  bool do_beamskip = self->do_beamskip_;
  double beam_skip_distance = self->beam_skip_distance_;
  double beam_skip_threshold = self->beam_skip_threshold_;

  // we only do beam skipping if the filter has converged
  if (do_beamskip && !set->converged) {
    do_beamskip = false;
  }

  // we need a count the no of particles for which the beam agreed with the map
  int *obs_count = new int[self->max_beams_]();

  // we also need a mask of which observations to integrate (to decide which
  // beams to integrate to all particles)
  bool *obs_mask = new bool[self->max_beams_]();

  int beam_ind = 0;

  // realloc indicates if we need to reallocate the temp data structure needed
  // to do beamskipping
  bool realloc = false;

  if (do_beamskip) {
    if (self->max_obs_ < self->max_beams_) {
      realloc = true;
    }

    if (self->max_samples_ < set->samples.size()) {
      realloc = true;
    }

    if (realloc) {
      self->reallocTempData(set->samples.size(), self->max_beams_);
      fprintf(stderr, "Reallocing temp weights %d - %d\n", self->max_samples_, self->max_obs_);
    }
  }

  // Compute the sample weights
  for (size_t j = 0; j < set->samples.size(); j++) {
    pf::SamplePtr sample = set->samples[j];
    pf::Pose pose = sample->pose;

    // Take account of the laser pose relative to the robot
    //pose = pf_vector_coord_add(self->laser_pose_, pose);
    pose = pf::Pose::CoordAdd(self->laser_pose_, pose);

    double log_p = 0;

    beam_ind = 0;

    for (size_t i = 0; i < data->ranges.size(); i += step, beam_ind++) {
      //obs_range = data->ranges[i][0];
      //obs_bearing = data->ranges[i][1];
      double obs_range;
      double obs_bearing = data->ranges[i].angle;

      // This model ignores max range readings
      if (obs_range >= data->range_max) {
        continue;
      }

      // Check for NaN
      if (obs_range != obs_range) {
        continue;
      }

      double pz = 0.0;

      pf::Pose hit;
      // Compute the endpoint of the beam
      hit.v[0] = pose.v[0] + obs_range * cos(pose.v[2] + obs_bearing);
      hit.v[1] = pose.v[1] + obs_range * sin(pose.v[2] + obs_bearing);

      // Convert to map grid coords.
      int mi, mj;
      mi = self->map_->toGridX(hit.v[0]);
      mj = self->map_->toGridY(hit.v[1]);

      // Part 1: Get distance from the hit to closest obstacle.
      // Off-map penalized as max distance

      //if (!MAP_VALID(self->map, mi, mj)) {
      if (self->map_->isValid(mi, mj)) {
        pz += self->z_hit_ * max_dist_prob;
      } else {
        double z = self->map_->cells[self->map_->toIndex(mi, mj)]->occ_dist;
        if (z < beam_skip_distance) {
          obs_count[beam_ind] += 1;
        }
        pz += self->z_hit_ * exp(-(z * z) / z_hit_denom);
      }

      // Gaussian model
      // NOTE: this should have a normalization of 1/(sqrt(2pi)*sigma)

      // Part 2: random measurements
      pz += self->z_rand_ * z_rand_mult;

      assert(pz <= 1.0);
      assert(pz >= 0.0);

      // TODO: outlier rejection for short readings

      if (!do_beamskip) {
        log_p += log(pz);
      } else {
        self->temp_obs[j][beam_ind] = pz;
      }
    }
    if (!do_beamskip) {
      sample->weight *= exp(log_p);
      total_weight += sample->weight;
    }
  }

  if (do_beamskip) {
    int skipped_beam_count = 0;
    for (beam_ind = 0; beam_ind < self->max_beams_; beam_ind++) {
      if ((obs_count[beam_ind] / static_cast<double>(set->samples.size())) >
          beam_skip_threshold) {
        obs_mask[beam_ind] = true;
      } else {
        obs_mask[beam_ind] = false;
        skipped_beam_count++;
      }
    }

    // we check if there is at least a critical number of beams that agreed with
    // the map
    // otherwise it probably indicates that the filter converged to a wrong
    // solution
    // if that's the case we integrate all the beams and hope the filter might
    // converge to
    // the right solution
    bool error = false;

    if (skipped_beam_count >= (beam_ind * self->beam_skip_error_threshold_)) {
      fprintf(stderr,
              "Over %f%% of the observations were not in the map - pf may have "
              "converged to wrong pose - integrating all observations\n",
              (100 * self->beam_skip_error_threshold_));
      error = true;
    }

    for (size_t j = 0; j < set->samples.size(); j++) {
      pf::SamplePtr sample = set->samples[j];
      pf::Pose pose = sample->pose;

      double log_p = 0;

      for (beam_ind = 0; beam_ind < self->max_beams_; beam_ind++) {
        if (error || obs_mask[beam_ind]) {
          log_p += log(self->temp_obs[j][beam_ind]);
        }
      }

      sample->weight *= exp(log_p);

      total_weight += sample->weight;
    }
  }

  delete[] obs_count;
  delete[] obs_mask;
  return (total_weight);
}

void AMCLLaser::reallocTempData(int new_max_samples, int new_max_obs) {
  if (temp_obs) {
    for (int k = 0; k < max_samples_; k++) {
      delete[] temp_obs[k];
    }
    delete[] temp_obs;
  }
  max_obs_ = new_max_obs;
  max_samples_ = fmax(max_samples_, new_max_samples);

  temp_obs = new double *[max_samples_]();
  for (int k = 0; k < max_samples_; k++) {
    temp_obs[k] = new double[max_obs_]();
  }
}
