/***
 * Desc: Simple particle filter for localization.
 * Author: QianHao Huang <qiao.helloworld@gmail.com>
 * Date: 8 Aug 2016
 */

#pragma once

#include "pf_vector.hpp"
#include "pf_kdtree.hpp"
//#include "../sensors/amcl_sensor.hpp"

#include <functional>
#include <memory>
#include <vector>


namespace amcl {
  struct AMCLSensorData;
  using AMCLSensorDataPtr = std::shared_ptr<AMCLSensorData>;
namespace map {
  struct Map;
  using MapPtr = std::shared_ptr<Map>;
} // namespace map
namespace pf {

// Forward declarations
struct Sample;
using SamplePtr = std::shared_ptr<Sample>;

struct SampleSet;
using SampleSetPtr = std::shared_ptr<SampleSet>;

struct Cluster;
using ClusterPtr = std::shared_ptr<Cluster>;

class ParticleFilter;
using ParticleFilterPtr = std::shared_ptr<ParticleFilter>;


// Function prototype for the initialization model; generates a sample pose from
// an appropriate distribution.
//typedef pf_vector_t (*pf_init_model_fn_t) (void *init_data);
using InitModelFunc = std::function<Pose(std::shared_ptr<map::Map> init_data)>;

// Function prototype for the action model; generates a sample pose from
// an appropriate distribution
//typedef void (*pf_action_model_fn_t) (void *action_data, struct _pf_sample_set_t* set);
using ActionModelFunc = std::function<void(void* action_data, SampleSet &set)>;

// Function prototype for the sensor model; determines the probability
// for the given set of sample poses.
//typedef double (*pf_sensor_model_fn_t) (void *sensor_data, struct _pf_sample_set_t* set);
using SensorModelFunc = std::function<double(amcl::AMCLSensorDataPtr sensor_data, pf::SampleSetPtr set)>;

// Information for a single sample
struct Sample
{
  // Pose represented by this sample
  Pose pose;

  // Weight for this pose
  double weight;
};// pf_sample_t;

// Information for a cluster of samples
struct Cluster
{
  // Number of samples
  int count;

  // Total weight of samples in this cluster
  double weight;

  // Cluster statistics
  Pose mean;
  Covariance cov;

  // Workspace
  double m[4], c[2][2];  
};// pf_cluster_t;

// Information for a set of samples
struct SampleSet //_pf_sample_set_t
{
  // The samples
  std::vector<SamplePtr> samples;

  // A kdtree encoding the histogram
  kdtree::Tree kdtree;

  // Clusters
  int cluster_count, cluster_max_count;
  std::vector<ClusterPtr> clusters;

  // Filter statistics
  Pose mean;
  Covariance cov;
  int converged; 
};// pf_sample_set_t;

// Information for an entire filter
class ParticleFilter // _pf_t
{
public:
  ParticleFilter(
          const size_t &min_samples, const size_t &max_samples,
          const double &alpha_slow, const double &alpha_fast,
          InitModelFunc random_pose_fn, amcl::map::MapPtr random_pose_data);

  inline size_t getMinSamples() { return min_samples_; }
  inline size_t getMaxSamples() { return max_samples_; }

  inline double getPopulationErr() { return pop_err_; }
  inline void setPopulationErr(const double &pop_err) { pop_err_ = pop_err; }

  inline double getPopulationZ() { return pop_z_; }
  inline void setPopulationZ(const double &pop_z) { pop_z_ = pop_z; }

  inline SampleSetPtr getCurrentSet() { return sets_[current_set_]; }

  // Initialize the filter using a guassian
  void initGuassian(const Pose &mean, const Covariance &cov);

  // Initailize the filter using some model
  void initModel(InitModelFunc init_fn, std::shared_ptr<map::Map> init_data);

  // Update the filter with some new action
  //void updateAction(ActionModelFunc, void* action_data);

  // Update the filter with some new sensor observation
  //void updateSensor(SensorModelFunc, void* sensor_data);

  // Resample the distribution
  void updateResample();

  // Compute the statistics for a particular cluster
  // Returns false if there is no such cluster.
  bool getClusterStats(const int &clabel, double& weight, Pose& mean, Covariance& cov);

  // Calculate if the particle filter has converged - and sets the converged
  // flag in the current set and the pf
  bool updateConverged();

  // Sets the current set and pf converged values to zero
  void resetConverge();

  void updateSampleWeights(const double &total_weight);

private:
  // Compute the required number of samples, given that there are k bins with samples in them.
  size_t resampleLimit(const size_t &k);

  // Re-compute the cluster statistics for a sample set
  void clusterStats(SampleSetPtr set);

private:
  // This min and max number of samples
  size_t min_samples_, max_samples_;

  // Population size parameters
  double pop_err_, pop_z_;

  // The sample sets.  We keep two sets and use [current_set]
  // to identify the active set.
  size_t current_set_;
  SampleSetPtr sets_[2];

  // Running averages, slow and fast, of likelihood
  double w_slow_, w_fast_;

  // Decay rates for running averages
  double alpha_slow_, alpha_fast_;

  // Function used to draw random pose samples
  InitModelFunc random_pose_fn_;
  //void *random_pose_data_;
  //std::unique_ptr<void> random_pose_data_;
  std::shared_ptr<map::Map> random_pose_data_;

  // distance threshold in each axis over which the pf is considered to not be converged
  double dist_threshold_;

  bool converged_;
};// pf_t;

} // namespace pf
} // namespace amcl
