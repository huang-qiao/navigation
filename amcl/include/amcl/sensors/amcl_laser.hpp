#pragma once

#include "amcl_sensor.hpp"
#include "../map/map.hpp"

namespace amcl {

struct AMCLLaserData;
using AMCLLaserDataPtr = std::shared_ptr<AMCLLaserData>;

class AMCLLaser;
using AMCLLaserPtr = std::shared_ptr<AMCLLaser>;

typedef enum
{
  LASER_MODEL_BEAM,
  LASER_MODEL_LIKELIHOOD_FIELD,
  LASER_MODEL_LIKELIHOOD_FIELD_PROB
} LaserModelType;

struct AMCLLaserBeam
{
  double angle;
  double range;
};

// Laser sensor data
struct AMCLLaserData : public AMCLSensorData
{
public:
  double angle_min;
  double angle_max;
  double angle_increment;

  double range_min;
  double range_max;

  std::vector<AMCLLaserBeam> ranges;
};


// Laseretric sensor model
class AMCLLaser : public AMCLSensor
{
public:
  // Default constructor
  AMCLLaser(const size_t &max_beams, const map::MapPtr map);

  virtual ~AMCLLaser();

  void setModelBeam(double z_hit,
                    double z_short,
                    double z_max,
                    double z_rand,
                    double sigma_hit,
                    double labda_short,
                    double chi_outlier);

  void setModelLikelihoodField(double z_hit,
                               double z_rand,
                               double sigma_hit,
                               double max_occ_dist);

  //a more probabilistically correct model - also with the option to do beam skipping
  void setModelLikelihoodFieldProb(double z_hit,
                                   double z_rand,
                                   double sigma_hit,
                                   double max_occ_dist,
                                   bool do_beamskip,
                                   double beam_skip_distance,
                                   double beam_skip_threshold,
                                   double beam_skip_error_threshold);

  // Update the filter based on the sensor model.  Returns true if the
  // filter has been updated.
  virtual bool updateSensor(pf::ParticleFilterPtr pf, AMCLLaserDataPtr data);

  // Set the laser's pose after construction
  inline void setLaserPose(pf::Pose &laser_pose) { laser_pose_ = laser_pose; }

private:
  // Determine the probability for the given pose
  static double BeamModel(AMCLLaserDataPtr data, pf::SampleSetPtr set);
  // Determine the probability for the given pose
  static double LikelihoodFieldModel(AMCLLaserDataPtr data, pf::SampleSetPtr set);

  // Determine the probability for the given pose - more probablistic model 
  static double LikelihoodFieldModelProb(AMCLLaserDataPtr data, pf::SampleSetPtr set);

  void reallocTempData(int max_samples, int max_obs);

  LaserModelType model_type_;

  // Current data timestamp
  double time_;

  // The laser map
  map::MapPtr map_;

  // Laser offset relative to robot
  pf::Pose laser_pose_;
  
  // Max beams to consider
  int max_beams_;

  // Beam skipping parameters (used by LikelihoodFieldModelProb model)
  bool do_beamskip_;
  double beam_skip_distance_;
  double beam_skip_threshold_;
  //threshold for the ratio of invalid beams - at which all beams are integrated to the likelihoods 
  //this would be an error condition 
  double beam_skip_error_threshold_;

  //temp data that is kept before observations are integrated to each particle (requried for beam skipping)
  int max_samples_;
  int max_obs_;
  double **temp_obs;

  // Laser model params
  //
  // Mixture params for the components of the model; must sum to 1
  double z_hit_;
  double z_short_;
  double z_max_;
  double z_rand_;
  //
  // Stddev of Gaussian model for laser hits.
  double sigma_hit_;
  // Decay rate of exponential model for short readings.
  double lambda_short_;
  // Threshold for outlier rejection (unused)
  double chi_outlier_;
};

} // namespace amcl

