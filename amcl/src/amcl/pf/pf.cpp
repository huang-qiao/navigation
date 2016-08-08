#include <cassert>
#include <cmath>
#include <cstdlib>
#include <ctime>

#include "pf.hpp"
#include "pf_pdf.hpp"
#include "pf_kdtree.hpp"

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

using namespace amcl::pf;

// Compute the required number of samples, given that there are k bins
// with samples in them.
// static int pf_resample_limit(pf_t *pf, int k);
size_t ParticleFilter::resampleLimit(const size_t &k) {
  if (k <= 1) return max_samples_;

  double a, b, c, x;
  size_t n;

  a = 1.0;
  b = 2.0 / (9.0 * ((double)k - 1.0));
  c = sqrt(2.0 / (9 * ((double)k - 1.0))) * pop_z_;
  x = a - b + c;

  n = (size_t)ceil((k - 1.0) / (2 * pop_err_) * x * x * x);

  if (n < min_samples_) return min_samples_;
  if (n > max_samples_) return max_samples_;

  return n;
}

// Re-compute the cluster statistics for a sample set
// static void pf_cluster_stats(pf_t *pf, pf_sample_set_t *set);
void ParticleFilter::clusterStats(SampleSetPtr set) {
  TRACE_FUNC_ENTER
  // int cidx;

  // Cluster the samples
  // pf_kdtree_cluster(set->kdtree);
  set->kdtree.cluster();

  // Initialize cluster stats
  set->cluster_count = 0;
  set->clusters.resize(set->cluster_max_count);
  for (size_t i = 0; i < set->cluster_max_count; i++) {
    set->clusters[i] = std::make_shared<Cluster>();
    ClusterPtr cluster = set->clusters[i];
    cluster->count = 0;
    cluster->weight = 0;
    cluster->mean = Pose();
    cluster->cov = Covariance();

    for (size_t i = 0; i < 4; i++) {
      cluster->m[i] = 0.0;
    }
    for (size_t i = 0; i < 2; i++) {
      for (size_t j = 0; j < 2; j++) {
        cluster->c[i][j] = 0.0;
      }
    }
  }
  // Initialize overall filter stats
  size_t count = 0;
  double weight = 0.0;
  double m[4], c[2][2];

  set->mean = Pose();
  set->cov = Covariance();

  for (size_t i = 0; i < 4; i++) {
    m[i] = 0.0;
  }
  for (size_t i = 0; i < 2; i++) {
    for (size_t j = 0; j < 2; j++) {
      c[i][j] = 0.0;
    }
  }
  // Compute cluster stats
  for (size_t i = 0; i < set->samples.size(); i++) {
    SamplePtr sample = set->samples[i];

    // printf("%d %f %f %f\n", i, sample->pose.v[0], sample->pose.v[1],
    // sample->pose.v[2]);

    // Get the cluster label for this sample
    // cidx = pf_kdtree_get_cluster(set->kdtree, sample->pose);
    size_t cidx = set->kdtree.getCluster(sample->pose);
    assert(cidx >= 0);

    if (cidx >= set->cluster_max_count) {
      continue;
    }
    if (cidx + 1 > set->cluster_count) {
      set->cluster_count = cidx + 1;
    }

    set->clusters[cidx]->count += 1;
    set->clusters[cidx]->weight += sample->weight;

    count += 1;
    weight += sample->weight;

    // Compute mean
    set->clusters[cidx]->m[0] += sample->weight * sample->pose.v[0];
    set->clusters[cidx]->m[1] += sample->weight * sample->pose.v[1];
    set->clusters[cidx]->m[2] += sample->weight * cos(sample->pose.v[2]);
    set->clusters[cidx]->m[3] += sample->weight * sin(sample->pose.v[2]);

    m[0] += sample->weight * sample->pose.v[0];
    m[1] += sample->weight * sample->pose.v[1];
    m[2] += sample->weight * cos(sample->pose.v[2]);
    m[3] += sample->weight * sin(sample->pose.v[2]);

    // Compute covariance in linear components
    for (size_t i = 0; i < 2; i++)
      for (size_t j = 0; j < 2; j++) {
        set->clusters[cidx]->c[i][j] +=
            sample->weight * sample->pose.v[i] * sample->pose.v[j];
        c[i][j] += sample->weight * sample->pose.v[i] * sample->pose.v[j];
      }
  }
  // Normalize
  for (size_t i = 0; i < set->cluster_count; i++) {
    // cluster = set->clusters + i;

    set->clusters[i]->mean.v[0] =
        set->clusters[i]->m[0] / set->clusters[i]->weight;
    set->clusters[i]->mean.v[1] =
        set->clusters[i]->m[1] / set->clusters[i]->weight;
    set->clusters[i]->mean.v[2] =
        atan2(set->clusters[i]->m[3], set->clusters[i]->m[2]);

    set->clusters[i]->cov = Covariance();

    // Covariance in linear components
    for (size_t j = 0; j < 2; j++) {
      for (size_t k = 0; k < 2; k++) {
        set->clusters[i]->cov.m[j][k] =
            set->clusters[i]->c[j][k] / set->clusters[i]->weight -
            set->clusters[i]->mean.v[j] * set->clusters[i]->mean.v[k];
      }
    }
    // Covariance in angular components; I think this is the correct
    // formula for circular statistics.
    set->clusters[i]->cov.m[2][2] =
        -2 * log(sqrt(set->clusters[i]->m[2] * set->clusters[i]->m[2] +
                      set->clusters[i]->m[3] * set->clusters[i]->m[3]));

    // printf("cluster %d %d %f (%f %f %f)\n", i, cluster->count,
    // cluster->weight,
    // cluster->mean.v[0], cluster->mean.v[1], cluster->mean.v[2]);
    // pf_matrix_fprintf(cluster->cov, stdout, "%e");
  }
  // Compute overall filter stats
  set->mean.v[0] = m[0] / weight;
  set->mean.v[1] = m[1] / weight;
  set->mean.v[2] = atan2(m[3], m[2]);
  // Covariance in linear components
  for (size_t j = 0; j < 2; j++) {
    for (size_t k = 0; k < 2; k++) {
      set->cov.m[j][k] = c[j][k] / weight - set->mean.v[j] * set->mean.v[k];
    }
  }
  // Covariance in angular components; I think this is the correct
  // formula for circular statistics.
  set->cov.m[2][2] = -2 * log(sqrt(m[2] * m[2] + m[3] * m[3]));
}

// Create a new filter
// pf_t *pf_alloc(int min_samples, int max_samples,
//               double alpha_slow, double alpha_fast,
//               pf_init_model_fn_t random_pose_fn, void *random_pose_data)
ParticleFilter::ParticleFilter(const size_t &min_samples,
                               const size_t &max_samples,
                               const double &alpha_slow,
                               const double &alpha_fast,
                               InitModelFunc random_pose_fn,
                               amcl::map::MapPtr random_pose_data)
    : sets_() {
  TRACE_FUNC_ENTER
  // int i, j;
  // pf_t *pf;
  // pf_sample_set_t *set;
  // pf_sample_t *sample;

  srand48(time(NULL));
  TRACE_FUNC

  // pf = calloc(1, sizeof(pf_t));

  random_pose_fn_ = random_pose_fn;
  random_pose_data_ = random_pose_data;

  min_samples_ = min_samples;
  max_samples_ = max_samples;

  // Control parameters for the population size calculation.
  // [err] the max error between the true distribution and the estimated
  // distribution.
  // [z] the upper standard normal quantile for (1 - p), where p is the
  // probability that the error on the estimated distrubition will be less than
  // [err].
  pop_err_ = 0.01;
  pop_z_ = 3;
  dist_threshold_ = 0.5;

  current_set_ = 0;
  TRACE_FUNC
  for (size_t j = 0; j < 2; j++) {
    // set = pf->sets + j;
    TRACE_FUNC
    // sets_[j]->sample_count = max_samples;
    // sets_[j]->samples = calloc(max_samples, sizeof(pf_sample_t));
    sets_[j] = std::make_shared<SampleSet>();
    sets_[j]->samples.resize(max_samples_);
    TRACE_FUNC
    for (size_t k = 0; k < sets_[j]->samples.size(); k++) {
      sets_[j]->samples[k] = std::make_shared<Sample>();
    }

    // for (i = 0; i < sets_[j]->sample_count; i++)
    for (size_t i = 0; i < sets_[j]->samples.size(); i++) {
      // sample = set->samples + i;
      sets_[j]->samples[i]->pose.v[0] = 0.0;
      sets_[j]->samples[i]->pose.v[1] = 0.0;
      sets_[j]->samples[i]->pose.v[2] = 0.0;
      sets_[j]->samples[i]->weight = 1.0 / max_samples_;
    }
    TRACE_FUNC
    // HACK: is 3 times max_samples enough?
    // sets_[j]->kdtree = pf_kdtree_alloc(3 * max_samples);
    sets_[j]->kdtree = *kdtree::Tree::CreateTree(3 * max_samples_);

    sets_[j]->cluster_count = 0;
    sets_[j]->cluster_max_count = max_samples;
    // sets_[j]->clusters = calloc(sets_[j]->cluster_max_count,
    // sizeof(Cluster));
    sets_[j]->clusters.resize(sets_[j]->cluster_max_count);

    sets_[j]->mean = Pose();
    sets_[j]->cov = Covariance();
  }
  TRACE_FUNC
  w_slow_ = 0.0;
  w_fast_ = 0.0;

  alpha_slow_ = alpha_slow;
  alpha_fast_ = alpha_fast;

  // set converged to 0
  // pf_init_converged(pf);
  converged_ = false;
  TRACE_FUNC_EXIT
}

/* Free an existing filter
void pf_free(pf_t *pf)
{
  int i;

  for (i = 0; i < 2; i++)
  {
    free(pf->sets[i].clusters);
    pf_kdtree_free(pf->sets[i].kdtree);
    free(pf->sets[i].samples);
  }
  free(pf);

  return;
}
*/

// Initialize the filter using a guassian
// void pf_init(pf_t *pf, pf_vector_t mean, pf_matrix_t cov)
void ParticleFilter::initGuassian(const Pose &mean, const Covariance &cov) {
  TRACE_FUNC_ENTER
  // Create the kd tree for adaptive sampling
  // pf_kdtree_clear(set->kdtree);
  sets_[current_set_]->kdtree.clear();
  TRACE_FUNC

  // set->sample_count = pf->max_samples;
  sets_[current_set_]->samples.resize(max_samples_);

  // pdf = pf_pdf_gaussian_alloc(mean, cov);
  std::cout << "CreatePdfGaussian: mean("
            << mean.v[0] << ","
            << mean.v[1] << ","
            << mean.v[2] << "), cov("
            << cov.m[0][0] << ","
            << cov.m[1][1] << ","
            << cov.m[2][2] << ")" << std::endl;
  PdfGaussianPtr pdf = PdfGaussian::CreatePdfGaussian(mean, cov);
  TRACE_FUNC

  // Compute the new sample poses
  for (size_t i = 0; i < sets_[current_set_]->samples.size(); i++) {
    // sample = set->samples + i;
    sets_[current_set_]->samples[i]->weight = 1.0 / max_samples_;
    sets_[current_set_]->samples[i]->pose =
        pdf->genGaussianSample();  // pf_pdf_gaussian_sample(pdf);

    // std::cout << "insert to kdtree: pose("
    //          << sets_[current_set_]->samples[i]->pose.v[0] << ","
    //          << sets_[current_set_]->samples[i]->pose.v[1] << ","
    //          << sets_[current_set_]->samples[i]->pose.v[2] << "), weight("
    //          << sets_[current_set_]->samples[i]->weight << ")" << std::endl;

    // Add sample to histogram
    // pf_kdtree_insert(set->kdtree, sample->pose, sample->weight);
    sets_[current_set_]->kdtree.insert(sets_[current_set_]->samples[i]->pose,
                                       sets_[current_set_]->samples[i]->weight);
  }

  w_slow_ = w_fast_ = 0.0;

  // pf_pdf_gaussian_free(pdf);
  pdf.reset();
  TRACE_FUNC

  // Re-compute cluster statistics
  // pf_cluster_stats(pf, set);
  clusterStats(sets_[current_set_]);
  TRACE_FUNC

  // set converged to 0
  // pf_init_converged(pf);
  converged_ = false;
  TRACE_FUNC_EXIT
}

// Initialize the filter using some model
// void pf_init_model(pf_t *pf, pf_init_model_fn_t init_fn, void *init_data)
void ParticleFilter::initModel(InitModelFunc init_fn,
                               std::shared_ptr<map::Map> init_data) {
  TRACE_FUNC_ENTER
  // int i;
  // pf_sample_set_t *set;
  // pf_sample_t *sample;

  SampleSetPtr set = sets_[current_set_];

  // Create the kd tree for adaptive sampling
  // pf_kdtree_clear(set->kdtree);
  set->kdtree.clear();

  // set->sample_count = pf->max_samples;
  set->samples.resize(max_samples_);

  // Compute the new sample poses
  for (size_t i = 0; i < set->samples.size(); i++) {
    SamplePtr sample = set->samples[i];
    sample->weight = 1.0 / max_samples_;
    // sample->pose = (*init_fn) (init_data);
    sample->pose = init_fn(init_data);

    // Add sample to histogram
    // pf_kdtree_insert(set->kdtree, sample->pose, sample->weight);
    set->kdtree.insert(sample->pose, sample->weight);
  }

  w_slow_ = w_fast_ = 0.0;

  // Re-compute cluster statistics
  // pf_cluster_stats(pf, set);
  clusterStats(set);

  // set converged to 0
  // pf_init_converged(pf);
  resetConverge();
  TRACE_FUNC_EXIT
  return;
}

// void pf_init_converged(pf_t *pf)
void ParticleFilter::resetConverge() {
  TRACE_FUNC_ENTER
  // pf_sample_set_t *set;
  SampleSetPtr set;
  set = sets_[current_set_];
  set->converged = 0;
  converged_ = 0;
  TRACE_FUNC_EXIT
}

// int pf_update_converged(pf_t *pf)
bool ParticleFilter::updateConverged() {
  TRACE_FUNC_ENTER
  // int i;
  // pf_sample_set_t *set;
  // pf_sample_t *sample;
  // double total;

  SampleSetPtr set = sets_[current_set_];
  double mean_x = 0, mean_y = 0;

  for (size_t i = 0; i < set->samples.size(); i++) {
    SamplePtr sample = set->samples[i];

    mean_x += sample->pose.v[0];
    mean_y += sample->pose.v[1];
  }
  mean_x /= (double)set->samples.size();
  mean_y /= (double)set->samples.size();

  for (size_t i = 0; i < set->samples.size(); i++) {
    SamplePtr sample = set->samples[i];
    if (fabs(sample->pose.v[0] - mean_x) > dist_threshold_ ||
        fabs(sample->pose.v[1] - mean_y) > dist_threshold_) {
      set->converged = false;
      converged_ = false;
      TRACE_FUNC_EXIT
      return false;
    }
  }
  set->converged = true;
  converged_ = true;
  TRACE_FUNC_EXIT
  return true;
}

/* Update the filter with some new action
void pf_update_action(pf_t *pf, pf_action_model_fn_t action_fn, void
*action_data)
{
  pf_sample_set_t *set;

  set = pf->sets + pf->current_set;

  (*action_fn) (action_data, set);

  return;
}
*/

void ParticleFilter::updateSampleWeights(const double &total_weight) {
  TRACE_FUNC_ENTER
  SampleSetPtr set = sets_[current_set_];

  if (total_weight > 0.0) {
    // Normalize weights
    double w_avg = 0.0;
    for (size_t i = 0; i < set->samples.size(); i++) {
      SamplePtr sample = set->samples[i];
      w_avg += sample->weight;
      sample->weight /= total_weight;
    }

    // Update running averages of likelihood of samples (Prob Rob p258)

    w_avg /= (double)set->samples.size();

    if (w_slow_ == 0.0) {
      w_slow_ = w_avg;
    } else {
      w_slow_ += alpha_slow_ * (w_avg - w_slow_);
    }

    if (w_fast_ == 0.0) {
      w_fast_ = w_avg;
    } else {
      w_fast_ += alpha_fast_ * (w_avg - w_fast_);
    }
    // printf("w_avg: %e slow: %e fast: %e\n",
    // w_avg, pf->w_slow, pf->w_fast);
  } else {
    // total_weight <= 0, Handle zero total
    for (size_t i = 0; i < set->samples.size(); i++) {
      SamplePtr sample = set->samples[i];
      sample->weight = 1.0 / set->samples.size();
    }
  }
  TRACE_FUNC_EXIT
}

/* Update the filter with some new sensor observation
#include <float.h>
void pf_update_sensor(pf_t *pf, pf_sensor_model_fn_t sensor_fn, void
*sensor_data)
{
  int i;
  pf_sample_set_t *set;
  pf_sample_t *sample;
  double total;

  set = pf->sets + pf->current_set;

  // Compute the sample weights
  total = (*sensor_fn) (sensor_data, set);

  if (total > 0.0)
  {
    // Normalize weights
    double w_avg=0.0;
    for (i = 0; i < set->sample_count; i++)
    {
      sample = set->samples + i;
      w_avg += sample->weight;
      sample->weight /= total;
    }
    // Update running averages of likelihood of samples (Prob Rob p258)
    w_avg /= set->sample_count;
    if(pf->w_slow == 0.0)
      pf->w_slow = w_avg;
    else
      pf->w_slow += pf->alpha_slow * (w_avg - pf->w_slow);
    if(pf->w_fast == 0.0)
      pf->w_fast = w_avg;
    else
      pf->w_fast += pf->alpha_fast * (w_avg - pf->w_fast);
    //printf("w_avg: %e slow: %e fast: %e\n",
           //w_avg, pf->w_slow, pf->w_fast);
  }
  else
  {
    // Handle zero total
    for (i = 0; i < set->sample_count; i++)
    {
      sample = set->samples + i;
      sample->weight = 1.0 / set->sample_count;
    }
  }

  return;
}
*/

// Resample the distribution
// void pf_update_resample(pf_t *pf)
void ParticleFilter::updateResample() {
  TRACE_FUNC_ENTER
  // int i;
  // double total;
  // pf_sample_set_t *set_a, *set_b;
  // pf_sample_t *sample_a, *sample_b;

  // double r,c,U;
  // int m;
  // double count_inv;
  // double* c;

  // double w_diff;

  SampleSetPtr set_a = sets_[current_set_];
  SampleSetPtr set_b = sets_[(current_set_ + 1) % 2];

  // Build up cumulative probability table for resampling.
  // TODO: Replace this with a more efficient procedure
  // (e.g.,
  // http://www.network-theory.co.uk/docs/gslref/GeneralDiscreteDistributions.html)
  // c = (double*)malloc(sizeof(double)*(set_a->sample_count+1));
  std::vector<double> c;
  c.resize(set_a->samples.size() + 1);
  c[0] = 0.0;
  for (size_t i = 0; i < set_a->samples.size(); i++) {
    c[i + 1] = c[i] + set_a->samples[i]->weight;
  }

  // Create the kd tree for adaptive sampling
  // pf_kdtree_clear(set_b->kdtree);
  set_b->kdtree.clear();

  // Draw samples from set a to create set b.
  double total = 0;
  // set_b->sample_count = 0;

  double w_diff = 1.0 - w_fast_ / w_slow_;
  if (w_diff < 0.0) {
    w_diff = 0.0;
  }
  // printf("w_diff: %9.6f\n", w_diff);

  // Can't (easily) combine low-variance sampler with KLD adaptive
  // sampling, so we'll take the more traditional route.
  /*
  // Low-variance resampler, taken from Probabilistic Robotics, p110
  count_inv = 1.0/set_a->sample_count;
  r = drand48() * count_inv;
  c = set_a->samples[0].weight;
  i = 0;
  m = 0;
  */
  SamplePtr sample_b;
  while (set_b->samples.size() < max_samples_) {
    // sample_b = set_b->samples + set_b->sample_count++;
    sample_b = std::make_shared<Sample>();

    if (drand48() < w_diff) {
      sample_b->pose = random_pose_fn_(random_pose_data_);
    } else {
      // Can't (easily) combine low-variance sampler with KLD adaptive
      // sampling, so we'll take the more traditional route.
      /*
      // Low-variance resampler, taken from Probabilistic Robotics, p110
      U = r + m * count_inv;
      while(U>c)
      {
        i++;
        // Handle wrap-around by resetting counters and picking a new random
        // number
        if(i >= set_a->sample_count)
        {
          r = drand48() * count_inv;
          c = set_a->samples[0].weight;
          i = 0;
          m = 0;
          U = r + m * count_inv;
          continue;
        }
        c += set_a->samples[i].weight;
      }
      m++;
      */

      // Naive discrete event sampler
      double r;
      r = drand48();
      size_t i;
      for (i = 0; i < set_a->samples.size(); i++) {
        if ((c[i] <= r) && (r < c[i + 1])) break;
      }
      assert(i < set_a->samples.size());

      SamplePtr sample_a = set_a->samples[i];

      assert(sample_a->weight > 0);

      // Add sample to list
      sample_b->pose = sample_a->pose;
    }

    sample_b->weight = 1.0;
    total += sample_b->weight;

    // Add sample to histogram
    // pf_kdtree_insert(set_b->kdtree, sample_b->pose, sample_b->weight);
    set_b->kdtree.insert(sample_b->pose, sample_b->weight);

    // See if we have enough samples yet
    if (set_b->samples.size() > resampleLimit(set_b->kdtree.leafCount())) {
      break;
    }
  }

  // Reset averages, to avoid spiraling off into complete randomness.
  if (w_diff > 0.0) {
    w_slow_ = w_fast_ = 0.0;
  }
  // fprintf(stderr, "\n\n");

  // Normalize weights
  for (size_t i = 0; i < set_b->samples.size(); i++) {
    sample_b = set_b->samples[i];
    sample_b->weight /= total;
  }

  // Re-compute cluster statistics
  // pf_cluster_stats(pf, set_b);
  clusterStats(set_b);

  // Use the newly created sample set
  current_set_ = (current_set_ + 1) % 2;

  // pf_update_converged(pf);
  updateConverged();

  // free(c);
  TRACE_FUNC_EXIT
  return;
}

/* Compute the CEP statistics (mean and variance).
void pf_get_cep_stats(pf_t *pf, pf_vector_t *mean, double *var)
{
  int i;
  double mn, mx, my, mrr;
  pf_sample_set_t *set;
  pf_sample_t *sample;

  set = pf->sets + pf->current_set;

  mn = 0.0;
  mx = 0.0;
  my = 0.0;
  mrr = 0.0;

  for (i = 0; i < set->sample_count; i++)
  {
    sample = set->samples + i;

    mn += sample->weight;
    mx += sample->weight * sample->pose.v[0];
    my += sample->weight * sample->pose.v[1];
    mrr += sample->weight * sample->pose.v[0] * sample->pose.v[0];
    mrr += sample->weight * sample->pose.v[1] * sample->pose.v[1];
  }

  mean->v[0] = mx / mn;
  mean->v[1] = my / mn;
  mean->v[2] = 0.0;

  *var = mrr / mn - (mx * mx / (mn * mn) + my * my / (mn * mn));

  return;
}*/

// Get the statistics for a particular cluster.
// int pf_get_cluster_stats(pf_t *pf, int clabel, double *weight, pf_vector_t
// *mean, pf_matrix_t *cov)
bool ParticleFilter::getClusterStats(const int &clabel, double &weight,
                                     Pose &mean, Covariance &cov) {
  TRACE_FUNC_ENTER
  SampleSetPtr set = sets_[current_set_];

  if (clabel >= set->cluster_count) {
    return false;
  }

  ClusterPtr cluster = set->clusters[clabel];

  weight = cluster->weight;
  mean = cluster->mean;
  cov = cluster->cov;
  TRACE_FUNC_EXIT
  return true;
}
