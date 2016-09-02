#include <cassert>
#include <cmath>
#include <cstdlib>
#include <ctime>

#include "pf.hpp"
#include "pf_kdtree.hpp"
#include "pf_pdf.hpp"

// Get the statistics for a particular cluster.
int SampleSet::getClusterStats(int clabel, double *weight, Pose *mean,
                                    Covariance *cov) {
  Cluster *cluster;


  if (clabel >= cluster_count)
    return 0;

  cluster = clusters + clabel;

  *weight = cluster->weight;
  *mean = cluster->mean;
  *cov = cluster->cov;

  return 1;
}

// Create a new filter
ParticleFilter::ParticleFilter(int min_samples, int max_samples,
                               double alpha_slow, double alpha_fast,
                               pf_init_model_fn_t random_pose_fn,
                               void *random_pose_data) {
  int i, j;
  // ParticleFilter *pf;
  SampleSet *set;
  Sample *sample;

  srand48(time(NULL));

  // pf = (ParticleFilter *)malloc(sizeof(ParticleFilter));

  random_pose_fn_ = random_pose_fn;
  random_pose_data_ = random_pose_data;

  min_samples_ = min_samples;
  max_samples_ = max_samples;

  // Control parameters for the population size calculation.  [err] is
  // the max error between the true distribution and the estimated
  // distribution.  [z] is the upper standard normal quantile for (1 -
  // p), where p is the probability that the error on the estimated
  // distrubition will be less than [err].
  pop_err_ = 0.01;
  pop_z_ = 3;
  dist_threshold_ = 0.5;

  current_set_ = 0;
  for (j = 0; j < 2; j++) {
    set = sets + j;

    set->sample_count = max_samples;
    set->samples = (Sample *)malloc(max_samples * sizeof(Sample));

    for (i = 0; i < set->sample_count; i++) {
      sample = set->samples + i;
      sample->pose.v[0] = 0.0;
      sample->pose.v[1] = 0.0;
      sample->pose.v[2] = 0.0;
      sample->weight = 1.0 / max_samples;
    }

    // HACK: is 3 times max_samples enough?
    // set->kdtree = pf_kdtree_alloc(3 * max_samples);
    set->kdtree = new KdTree(3 * max_samples);

    set->cluster_count = 0;
    set->cluster_max_count = max_samples;
    set->clusters = (Cluster *)malloc(set->cluster_max_count * sizeof(Cluster));

    set->mean = Pose(0.0, 0.0, 0.0); // set->mean = pf_vector_zero();
    set->cov = Covariance();         // pf_matrix_zero();
  }

  w_slow_ = 0.0;
  w_fast_ = 0.0;

  alpha_slow_ = alpha_slow;
  alpha_fast_ = alpha_fast;

  // set converged to 0
  initConverged();
}

ParticleFilter::~ParticleFilter() {
  for (size_t i = 0; i < 2; i++) {
    free(sets[i].clusters);
    delete sets[i].kdtree;
    free(sets[i].samples);
  }
}

/* Free an existing filter
void pf_delete(ParticleFilter *pf) {
  int i;

  for (i = 0; i < 2; i++) {
    free(pf->sets[i].clusters);
    // pf_kdtree_delete(pf->sets[i].kdtree);
    delete pf->sets[i].kdtree;
    free(pf->sets[i].samples);
  }
  free(pf);

  return;
}
*/
// Initialize the filter using a guassian
void ParticleFilter::init(Pose mean, Covariance cov) {
  int i;
  SampleSet *set;
  Sample *sample;
  PdfGaussianPtr pdf;

  set = sets + current_set_;

  // Create the kd tree for adaptive sampling
  // pf_kdtree_clear(set->kdtree);
  set->kdtree->clear();

  set->sample_count = max_samples_;

  pdf = PdfGaussian::CreatePdf(mean, cov);

  // Compute the new sample poses
  for (i = 0; i < set->sample_count; i++) {
    sample = set->samples + i;
    sample->weight = 1.0 / max_samples_;
    sample->pose = pdf->sample(); // sample->pose = PdfGaussian::sample(pdf);

    // Add sample to histogram
    // pf_kdtree_insert(set->kdtree, sample->pose, sample->weight);
    set->kdtree->insert(sample->pose, sample->weight);
  }

  w_slow_ = w_fast_ = 0.0;

  pdf.reset(); // pf_pdf_gaussian_free(pdf);

  // Re-compute cluster statistics
  clusterStats(set);

  // set converged to 0
  initConverged();
}

// Initialize the filter using some model
void ParticleFilter::initModel(pf_init_model_fn_t init_fn, void *init_data) {
  int i;
  SampleSet *set;
  Sample *sample;

  set = sets + current_set_;

  // Create the kd tree for adaptive sampling
  // pf_kdtree_clear(set->kdtree);
  set->kdtree->clear();

  set->sample_count = max_samples_;

  // Compute the new sample poses
  for (i = 0; i < set->sample_count; i++) {
    sample = set->samples + i;
    sample->weight = 1.0 / max_samples_;
    sample->pose = (*init_fn)(init_data);

    // Add sample to histogram
    // pf_kdtree_insert(set->kdtree, sample->pose, sample->weight);
    set->kdtree->insert(sample->pose, sample->weight);
  }

  w_slow_ = w_fast_ = 0.0;

  // Re-compute cluster statistics
  clusterStats(set);

  // set converged to 0
  initConverged();

  return;
}

void ParticleFilter::initConverged() {
  SampleSet *set;
  set = sets + current_set_;
  set->converged = 0;
  converged_ = 0;
}

int ParticleFilter::updateConverged() {
  int i;
  SampleSet *set;
  Sample *sample;
  double total;

  set = sets + current_set_;
  double mean_x = 0, mean_y = 0;

  for (i = 0; i < set->sample_count; i++) {
    sample = set->samples + i;

    mean_x += sample->pose.v[0];
    mean_y += sample->pose.v[1];
  }
  mean_x /= set->sample_count;
  mean_y /= set->sample_count;

  for (i = 0; i < set->sample_count; i++) {
    sample = set->samples + i;
    if (fabs(sample->pose.v[0] - mean_x) > dist_threshold_ ||
        fabs(sample->pose.v[1] - mean_y) > dist_threshold_) {
      set->converged = 0;
      converged_ = 0;
      return 0;
    }
  }
  set->converged = 1;
  converged_ = 1;
  return 1;
}

void ParticleFilter::normalizeWeights(const double &total_weight) {
  SampleSet *set;
  Sample *sample;
  double total = total_weight;
  int i;

  set = sets + current_set_;

  if (total > 0.0) {
    // Normalize weights
    double w_avg = 0.0;
    for (i = 0; i < set->sample_count; i++) {
      sample = set->samples + i;
      w_avg += sample->weight;
      sample->weight /= total;
    }
    // Update running averages of likelihood of samples (Prob Rob p258)
    w_avg /= set->sample_count;
    if (w_slow_ == 0.0)
      w_slow_ = w_avg;
    else
      w_slow_ += alpha_slow_ * (w_avg - w_slow_);
    if (w_fast_ == 0.0)
      w_fast_ = w_avg;
    else
      w_fast_ += alpha_fast_ * (w_avg - w_fast_);
    // printf("w_avg: %e slow: %e fast: %e\n",
    // w_avg, w_slow, w_fast);
  } else {
    // Handle zero total
    for (i = 0; i < set->sample_count; i++) {
      sample = set->samples + i;
      sample->weight = 1.0 / set->sample_count;
    }
  }
}

// Resample the distribution
void ParticleFilter::updateResample() {
  int i;
  double total;
  SampleSet *set_a, *set_b;
  Sample *sample_a, *sample_b;

  // double r,c,U;
  // int m;
  // double count_inv;
  double *c;

  double w_diff;

  set_a = sets + current_set_;
  set_b = sets + (current_set_ + 1) % 2;

  // Build up cumulative probability table for resampling.
  // TODO: Replace this with a more efficient procedure
  // (e.g.,
  // http://www.network-theory.co.uk/docs/gslref/GeneralDiscreteDistributions.html)
  c = (double *)malloc(sizeof(double) * (set_a->sample_count + 1));
  c[0] = 0.0;
  for (i = 0; i < set_a->sample_count; i++)
    c[i + 1] = c[i] + set_a->samples[i].weight;

  // Create the kd tree for adaptive sampling
  // pf_kdtree_clear(set_b->kdtree);
  set_b->kdtree->clear();

  // Draw samples from set a to create set b.
  total = 0;
  set_b->sample_count = 0;

  w_diff = 1.0 - w_fast_ / w_slow_;
  if (w_diff < 0.0)
    w_diff = 0.0;
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
  while (set_b->sample_count < max_samples_) {
    sample_b = set_b->samples + set_b->sample_count++;

    if (drand48() < w_diff)
      sample_b->pose = (random_pose_fn_)(random_pose_data_);
    else {
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
      for (i = 0; i < set_a->sample_count; i++) {
        if ((c[i] <= r) && (r < c[i + 1]))
          break;
      }
      assert(i < set_a->sample_count);

      sample_a = set_a->samples + i;

      assert(sample_a->weight > 0);

      // Add sample to list
      sample_b->pose = sample_a->pose;
    }

    sample_b->weight = 1.0;
    total += sample_b->weight;

    // Add sample to histogram
    // pf_kdtree_insert(set_b->kdtree, sample_b->pose, sample_b->weight);
    set_b->kdtree->insert(sample_b->pose, sample_b->weight);

    // See if we have enough samples yet
    if (set_b->sample_count > resampleLimit(set_b->kdtree->leaf_count))
      break;
  }

  // Reset averages, to avoid spiraling off into complete randomness.
  if (w_diff > 0.0)
    w_slow_ = w_fast_ = 0.0;

  // fprintf(stderr, "\n\n");

  // Normalize weights
  for (i = 0; i < set_b->sample_count; i++) {
    sample_b = set_b->samples + i;
    sample_b->weight /= total;
  }

  // Re-compute cluster statistics
  clusterStats(set_b);

  // Use the newly created sample set
  current_set_ = (current_set_ + 1) % 2;

  updateConverged();

  free(c);
  return;
}

// Compute the required number of samples, given that there are k bins
// with samples in them.  This is taken directly from Fox et al.
int ParticleFilter::resampleLimit(int k) {
  double a, b, c, x;
  int n;

  if (k <= 1)
    return max_samples_;

  a = 1;
  b = 2 / (9 * ((double)k - 1));
  c = sqrt(2 / (9 * ((double)k - 1))) * pop_z_;
  x = a - b + c;

  n = (int)ceil((k - 1) / (2 * pop_err_) * x * x * x);

  if (n < min_samples_)
    return min_samples_;
  if (n > max_samples_)
    return max_samples_;

  return n;
}

// Re-compute the cluster statistics for a sample set
void ParticleFilter::clusterStats(SampleSet *set) {
  int i, j, k, cidx;
  Sample *sample;
  Cluster *cluster;

  // Workspace
  double m[4], c[2][2];
  size_t count;
  double weight;

  // Cluster the samples
  // pf_kdtree_cluster(set->kdtree);
  set->kdtree->cluster();

  // Initialize cluster stats
  set->cluster_count = 0;

  for (i = 0; i < set->cluster_max_count; i++) {
    cluster = set->clusters + i;
    cluster->count = 0;
    cluster->weight = 0;
    cluster->mean = Pose(0.0, 0.0, 0.0); // pf_vector_zero();
    cluster->cov = Covariance();         // pf_matrix_zero();

    for (j = 0; j < 4; j++)
      cluster->m[j] = 0.0;
    for (j = 0; j < 2; j++)
      for (k = 0; k < 2; k++)
        cluster->c[j][k] = 0.0;
  }

  // Initialize overall filter stats
  count = 0;
  weight = 0.0;
  set->mean = Pose(0.0, 0.0, 0.0); // pf_vector_zero();
  set->cov = Covariance();         // pf_matrix_zero();
  for (j = 0; j < 4; j++)
    m[j] = 0.0;
  for (j = 0; j < 2; j++)
    for (k = 0; k < 2; k++)
      c[j][k] = 0.0;

  // Compute cluster stats
  for (i = 0; i < set->sample_count; i++) {
    sample = set->samples + i;

    // printf("%d %f %f %f\n", i, sample->pose.v[0], sample->pose.v[1],
    // sample->pose.v[2]);

    // Get the cluster label for this sample
    // cidx = pf_kdtree_get_cluster(set->kdtree, sample->pose);
    cidx = set->kdtree->getCluster(sample->pose);
    assert(cidx >= 0);
    if (cidx >= set->cluster_max_count)
      continue;
    if (cidx + 1 > set->cluster_count)
      set->cluster_count = cidx + 1;

    cluster = set->clusters + cidx;

    cluster->count += 1;
    cluster->weight += sample->weight;

    count += 1;
    weight += sample->weight;

    // Compute mean
    cluster->m[0] += sample->weight * sample->pose.v[0];
    cluster->m[1] += sample->weight * sample->pose.v[1];
    cluster->m[2] += sample->weight * cos(sample->pose.v[2]);
    cluster->m[3] += sample->weight * sin(sample->pose.v[2]);

    m[0] += sample->weight * sample->pose.v[0];
    m[1] += sample->weight * sample->pose.v[1];
    m[2] += sample->weight * cos(sample->pose.v[2]);
    m[3] += sample->weight * sin(sample->pose.v[2]);

    // Compute covariance in linear components
    for (j = 0; j < 2; j++)
      for (k = 0; k < 2; k++) {
        cluster->c[j][k] +=
            sample->weight * sample->pose.v[j] * sample->pose.v[k];
        c[j][k] += sample->weight * sample->pose.v[j] * sample->pose.v[k];
      }
  }

  // Normalize
  for (i = 0; i < set->cluster_count; i++) {
    cluster = set->clusters + i;

    cluster->mean.v[0] = cluster->m[0] / cluster->weight;
    cluster->mean.v[1] = cluster->m[1] / cluster->weight;
    cluster->mean.v[2] = atan2(cluster->m[3], cluster->m[2]);

    cluster->cov = Covariance(); // pf_matrix_zero();

    // Covariance in linear components
    for (j = 0; j < 2; j++)
      for (k = 0; k < 2; k++)
        cluster->cov.m[j][k] = cluster->c[j][k] / cluster->weight -
                               cluster->mean.v[j] * cluster->mean.v[k];

    // Covariance in angular components; I think this is the correct
    // formula for circular statistics.
    cluster->cov.m[2][2] = -2 * log(sqrt(cluster->m[2] * cluster->m[2] +
                                         cluster->m[3] * cluster->m[3]));

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
  for (j = 0; j < 2; j++)
    for (k = 0; k < 2; k++)
      set->cov.m[j][k] = c[j][k] / weight - set->mean.v[j] * set->mean.v[k];

  // Covariance in angular components; I think this is the correct
  // formula for circular statistics.
  set->cov.m[2][2] = -2 * log(sqrt(m[2] * m[2] + m[3] * m[3]));

  return;
}

// Compute the CEP statistics (mean and variance).
void ParticleFilter::getCepStats(Pose *mean, double *var) {
  int i;
  double mn, mx, my, mrr;
  SampleSet *set;
  Sample *sample;

  set = sets + current_set_;

  mn = 0.0;
  mx = 0.0;
  my = 0.0;
  mrr = 0.0;

  for (i = 0; i < set->sample_count; i++) {
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
}

// Get the statistics for a particular cluster.
int ParticleFilter::getClusterStats(int clabel, double *weight, Pose *mean,
                                    Covariance *cov) {
  SampleSet *set;
  Cluster *cluster;

  set = sets + current_set_;

  /*
  if (clabel >= set->cluster_count)
    return 0;
  cluster = set->clusters + clabel;

  *weight = cluster->weight;
  *mean = cluster->mean;
  *cov = cluster->cov;
  */

  int ret = set->getClusterStats(clabel, weight, mean, cov);

  return ret;
}
