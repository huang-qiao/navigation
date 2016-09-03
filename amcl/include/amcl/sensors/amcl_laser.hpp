#pragma once

#include "amcl_sensor.hpp"
#include "../map/map.hpp"

namespace amcl
{

enum class LaserModel
{
  BEAM,
  LIKELIHOOD_FIELD,
  LIKELIHOOD_FIELD_PROB
};


// Laser sensor data
struct AMCLLaserData : public AMCLSensorData
{
  //AMCLLaserData () {ranges=NULL;};
  //virtual ~AMCLLaserData() {delete [] ranges;};

  // Laser range data (range, bearing tuples)
  //int range_count;
  double angle_min, angle_max, angle_increment;
  double range_min, range_max;
  //double (*ranges)[2];
  std::vector<double> ranges;
};

using AMCLLaserDataPtr = std::shared_ptr<AMCLLaserData>;

// Laseretric sensor model
class AMCLLaser : public AMCLSensor
{
public:
  // Default constructor
  AMCLLaser(size_t max_beams, MapPtr map);

  virtual ~AMCLLaser();

  void SetModelBeam(double z_hit,
                            double z_short,
                            double z_max,
                            double z_rand,
                            double sigma_hit,
                            double labda_short,
                            double chi_outlier);

  void SetModelLikelihoodField(double z_hit,
                                       double z_rand,
                                       double sigma_hit,
                                       double max_occ_dist);

  //a more probabilistically correct model - also with the option to do beam skipping
  void SetModelLikelihoodFieldProb(double z_hit,
                                           double z_rand,
                                           double sigma_hit,
                                           double max_occ_dist,
                                           bool do_beamskip,
                                           double beam_skip_distance,
                                           double beam_skip_threshold,
                                           double beam_skip_error_threshold);

  // Update the filter based on the sensor model.  Returns true if the
  // filter has been updated.
  virtual bool UpdateSensor(ParticleFilter *pf, AMCLSensorDataPtr data);
  //virtual bool UpdateSensor(pf_t *pf, AMCLLaserDataPtr data);

  // Set the laser's pose after construction
  inline void SetLaserPose(Pose& laser_pose)
          {this->laser_pose = laser_pose;}

  // The laser's frame id
  inline void SetLaserFrameId(std::string laser_frame_id) { laser_frame_id_ = laser_frame_id; }
  inline std::string GetLaserFrameId() { return laser_frame_id_; }

private:
  // Determine the probability for the given pose
  static double BeamModel(AMCLLaserDataPtr data,
                                   SampleSet* set);
  // Determine the probability for the given pose
  static double LikelihoodFieldModel(AMCLLaserDataPtr data,
                                              SampleSet* set);

  // Determine the probability for the given pose - more probablistic model
  static double LikelihoodFieldModelProb(AMCLLaserDataPtr data,
                                             SampleSet* set);

  void reallocTempData(int max_samples, int max_obs);

  LaserModel model_type;

  // Current data timestamp
  double time;

  // The laser map
  MapPtr map;

  // Laser offset relative to robot
  Pose laser_pose;

  // Laser id
  std::string laser_frame_id_;

  // Max beams to consider
  int max_beams;

  // Beam skipping parameters (used by LikelihoodFieldModelProb model)
  bool do_beamskip;
  double beam_skip_distance;
  double beam_skip_threshold;
  //threshold for the ratio of invalid beams - at which all beams are integrated to the likelihoods
  //this would be an error condition
  double beam_skip_error_threshold;

  //temp data that is kept before observations are integrated to each particle (requried for beam skipping)
  int max_samples;
  int max_obs;
  double **temp_obs;

  // Laser model params
  //
  // Mixture params for the components of the model; must sum to 1
  double z_hit;
  double z_short;
  double z_max;
  double z_rand;
  //
  // Stddev of Gaussian model for laser hits.
  double sigma_hit;
  // Decay rate of exponential model for short readings.
  double lambda_short;
  // Threshold for outlier rejection (unused)
  double chi_outlier;
};

using AMCLLaserPtr = std::shared_ptr<AMCLLaser>;

}
