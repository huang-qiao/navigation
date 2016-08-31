#pragma once

#include "pf_kdtree.hpp"
#include "pf_vector.hpp"

// Forward declarations
struct ParticleFilter;
struct SampleSet;

// Function prototype for the initialization model; generates a sample pose from
// an appropriate distribution.
typedef Pose (*pf_init_model_fn_t)(void *init_data);

// Information for a single sample
struct Sample {
  // Pose represented by this sample
  Pose pose;

  // Weight for this pose
  double weight;
};

// Information for a cluster of samples
struct Cluster {
  // Number of samples
  int count;

  // Total weight of samples in this cluster
  double weight;

  // Cluster statistics
  Pose mean;
  Covariance cov;

  // Workspace
  double m[4], c[2][2];
};

// Information for a set of samples
struct SampleSet {
  // The samples
  int sample_count;
  Sample *samples;

  // A kdtree encoding the histogram
  KdTree *kdtree;

  // Clusters
  int cluster_count, cluster_max_count;
  Cluster *clusters;

  // Filter statistics
  Pose mean;
  Covariance cov;
  int converged;
};

// Information for an entire filter
struct ParticleFilter {
  // This min and max number of samples
  int min_samples_, max_samples_;

  // Population size parameters
  double pop_err_, pop_z_;

  // The sample sets.  We keep two sets and use [current_set]
  // to identify the active set.
  int current_set_;
  SampleSet sets[2];

  // Running averages, slow and fast, of likelihood
  double w_slow_, w_fast_;

  // Decay rates for running averages
  double alpha_slow_, alpha_fast_;

  // Function used to draw random pose samples
  pf_init_model_fn_t random_pose_fn_;
  void *random_pose_data_;

  double
      dist_threshold_; // distance threshold in each axis over which the pf is
                       // considered to not be converged
  int converged_;

  ParticleFilter(int min_samples_, int max_samples_, double alpha_slow_,
                 double alpha_fast_, pf_init_model_fn_t random_pose_fn_,
                 void *random_pose_data_);
  virtual ~ParticleFilter();

  // Initialize the filter using a guassian
  void init(Pose mean, Covariance cov);

  // Initialize the filter using some model
  void initModel(pf_init_model_fn_t init_fn, void *init_data);

  // normalize the weights of all particles
  void normalizeWeights(const double &total_weight);

  // Resample the distribution
  void updateResample();

  // Compute the CEP statistics (mean and variance).
  void getCepStats(Pose *mean, double *var);

  // Compute the statistics for a particular cluster.  Returns 0 if
  // there is no such cluster.
  int getClusterStats(int cluster, double *weight, Pose *mean, Covariance *cov);

  // sets the current set and pf converged values to zero
  void initConverged();

  // calculate if the particle filter has converged -
  // and sets the converged flag in the current set and the pf
  int updateConverged();

private:
  // Compute the required number of samples, given that there are k bins
  // with samples in them.
  int resampleLimit(int k);

  // Re-compute the cluster statistics for a sample set
  void clusterStats(SampleSet *set);

#ifdef INCLUDE_RTKGUI
public:
  // Display the sample set
  void pf_draw_samples(ParticleFilter *pf, struct _rtk_fig_t *fig,
                       int max_samples);

  // Draw the histogram (kdtree)
  void pf_draw_hist(ParticleFilter *pf, struct _rtk_fig_t *fig);

  // Draw the CEP statistics
  void pf_draw_cep_stats(ParticleFilter *pf, struct _rtk_fig_t *fig);

  // Draw the cluster statistics
  void pf_draw_cluster_stats(ParticleFilter *pf, struct _rtk_fig_t *fig);
#endif
};
