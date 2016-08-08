#include <algorithm>

#include <sys/types.h>  // required by Darwin
#include <cmath>

#include "amcl_odom.hpp"

using namespace amcl;

static double normalize(double z) { return atan2(sin(z), cos(z)); }

static double angle_diff(double a, double b) {
  double d1, d2;
  a = normalize(a);
  b = normalize(b);
  d1 = a - b;
  d2 = 2 * M_PI - fabs(d1);
  if (d1 > 0) d2 *= -1.0;
  if (fabs(d1) < fabs(d2))
    return (d1);
  else
    return (d2);
}

////////////////////////////////////////////////////////////////////////////////
// Default constructor
AMCLOdom::AMCLOdom() : AMCLSensor() { time_ = 0.0; }

void AMCLOdom::setModelDiff(double alpha1, double alpha2, double alpha3,
                            double alpha4) {
  model_type_ = ODOM_MODEL_DIFF;
  alpha1_ = alpha1;
  alpha2_ = alpha2;
  alpha3_ = alpha3;
  alpha4_ = alpha4;
}

void AMCLOdom::setModelOmni(double alpha1, double alpha2, double alpha3,
                            double alpha4, double alpha5) {
  model_type_ = ODOM_MODEL_OMNI;
  alpha1_ = alpha1;
  alpha2_ = alpha2;
  alpha3_ = alpha3;
  alpha4_ = alpha4;
  alpha5_ = alpha5;
}

void AMCLOdom::setModel(OdomModelType type, double alpha1, double alpha2,
                        double alpha3, double alpha4, double alpha5) {
  model_type_ = type;
  alpha1_ = alpha1;
  alpha2_ = alpha2;
  alpha3_ = alpha3;
  alpha4_ = alpha4;
  alpha5_ = alpha5;
}

////////////////////////////////////////////////////////////////////////////////
// Apply the action model
bool AMCLOdom::updateAction(pf::ParticleFilterPtr pf, AMCLSensorDataPtr data) {
  // AMCLOdomData *ndata;
  // ndata = (AMCLOdomData*) data;
  AMCLOdomDataPtr ndata = std::dynamic_pointer_cast<AMCLOdomData>(data);

  // Compute the new sample poses
  // pf::SampleSetPtr set;

  // set = pf->sets + pf->current_set;
  pf::SampleSetPtr set = pf->getCurrentSet();
  // pf_vector_t old_pose = pf_vector_sub(ndata->pose, ndata->delta);
  pf::Pose old_pose = ndata->pose - ndata->delta;

  switch (model_type_) {
    case ODOM_MODEL_OMNI: {
      double delta_trans, delta_rot, delta_bearing;
      double delta_trans_hat, delta_rot_hat, delta_strafe_hat;
      delta_trans = sqrt(ndata->delta.v[0] * ndata->delta.v[0] +
                         ndata->delta.v[1] * ndata->delta.v[1]);
      delta_rot = ndata->delta.v[2];

      // Precompute a couple of things
      double trans_hat_stddev = (alpha3_ * (delta_trans * delta_trans) +
                                 alpha1_ * (delta_rot * delta_rot));
      double rot_hat_stddev = (alpha4_ * (delta_rot * delta_rot) +
                               alpha2_ * (delta_trans * delta_trans));
      double strafe_hat_stddev = (alpha1_ * (delta_rot * delta_rot) +
                                  alpha5_ * (delta_trans * delta_trans));

      for (size_t i = 0; i < set->samples.size(); i++) {
        pf::SamplePtr sample = set->samples[i];

        delta_bearing = angle_diff(atan2(ndata->delta.v[1], ndata->delta.v[0]),
                                   old_pose.v[2]) +
                        sample->pose.v[2];
        double cs_bearing = cos(delta_bearing);
        double sn_bearing = sin(delta_bearing);

        // Sample pose differences
        delta_trans_hat = delta_trans + pf::RandomGaussian(trans_hat_stddev);
        delta_rot_hat = delta_rot + pf::RandomGaussian(rot_hat_stddev);
        delta_strafe_hat = 0 + pf::RandomGaussian(strafe_hat_stddev);
        // Apply sampled update to particle pose
        sample->pose.v[0] +=
            (delta_trans_hat * cs_bearing + delta_strafe_hat * sn_bearing);
        sample->pose.v[1] +=
            (delta_trans_hat * sn_bearing - delta_strafe_hat * cs_bearing);
        sample->pose.v[2] += delta_rot_hat;
      }
    } break;
    case ODOM_MODEL_DIFF: {
      // Implement sample_motion_odometry (Prob Rob p 136)
      double delta_rot1, delta_trans, delta_rot2;
      double delta_rot1_hat, delta_trans_hat, delta_rot2_hat;
      double delta_rot1_noise, delta_rot2_noise;

      // Avoid computing a bearing from two poses that are extremely near each
      // other (happens on in-place rotation).
      if (sqrt(ndata->delta.v[1] * ndata->delta.v[1] +
               ndata->delta.v[0] * ndata->delta.v[0]) < 0.01)
        delta_rot1 = 0.0;
      else
        delta_rot1 = angle_diff(atan2(ndata->delta.v[1], ndata->delta.v[0]),
                                old_pose.v[2]);
      delta_trans = sqrt(ndata->delta.v[0] * ndata->delta.v[0] +
                         ndata->delta.v[1] * ndata->delta.v[1]);
      delta_rot2 = angle_diff(ndata->delta.v[2], delta_rot1);

      // We want to treat backward and forward motion symmetrically for the
      // noise model to be applied below.  The standard model seems to assume
      // forward motion.
      delta_rot1_noise = std::min(fabs(angle_diff(delta_rot1, 0.0)),
                                  fabs(angle_diff(delta_rot1, M_PI)));
      delta_rot2_noise = std::min(fabs(angle_diff(delta_rot2, 0.0)),
                                  fabs(angle_diff(delta_rot2, M_PI)));

      for (int i = 0; i < set->samples.size(); i++) {
        pf::SamplePtr sample = set->samples[i];

        // Sample pose differences
        delta_rot1_hat = angle_diff(
            delta_rot1,
            pf::RandomGaussian(alpha1_ * delta_rot1_noise * delta_rot1_noise +
                               alpha2_ * delta_trans * delta_trans));
        delta_trans_hat =
            delta_trans -
            pf::RandomGaussian(alpha3_ * delta_trans * delta_trans +
                               alpha4_ * delta_rot1_noise * delta_rot1_noise +
                               alpha4_ * delta_rot2_noise * delta_rot2_noise);
        delta_rot2_hat = angle_diff(
            delta_rot2,
            pf::RandomGaussian(alpha1_ * delta_rot2_noise * delta_rot2_noise +
                               alpha2_ * delta_trans * delta_trans));

        // Apply sampled update to particle pose
        sample->pose.v[0] +=
            delta_trans_hat * cos(sample->pose.v[2] + delta_rot1_hat);
        sample->pose.v[1] +=
            delta_trans_hat * sin(sample->pose.v[2] + delta_rot1_hat);
        sample->pose.v[2] += delta_rot1_hat + delta_rot2_hat;
      }
    } break;
    case ODOM_MODEL_OMNI_CORRECTED: {
      double delta_trans, delta_rot, delta_bearing;
      double delta_trans_hat, delta_rot_hat, delta_strafe_hat;

      delta_trans = sqrt(ndata->delta.v[0] * ndata->delta.v[0] +
                         ndata->delta.v[1] * ndata->delta.v[1]);
      delta_rot = ndata->delta.v[2];

      // Precompute a couple of things
      double trans_hat_stddev = sqrt(alpha3_ * (delta_trans * delta_trans) +
                                     alpha1_ * (delta_rot * delta_rot));
      double rot_hat_stddev = sqrt(alpha4_ * (delta_rot * delta_rot) +
                                   alpha2_ * (delta_trans * delta_trans));
      double strafe_hat_stddev = sqrt(alpha1_ * (delta_rot * delta_rot) +
                                      alpha5_ * (delta_trans * delta_trans));

      for (int i = 0; i < set->samples.size(); i++) {
        pf::SamplePtr sample = set->samples[i];

        delta_bearing = angle_diff(atan2(ndata->delta.v[1], ndata->delta.v[0]),
                                   old_pose.v[2]) +
                        sample->pose.v[2];
        double cs_bearing = cos(delta_bearing);
        double sn_bearing = sin(delta_bearing);

        // Sample pose differences
        delta_trans_hat = delta_trans + pf::RandomGaussian(trans_hat_stddev);
        delta_rot_hat = delta_rot + pf::RandomGaussian(rot_hat_stddev);
        delta_strafe_hat = 0 + pf::RandomGaussian(strafe_hat_stddev);
        // Apply sampled update to particle pose
        sample->pose.v[0] +=
            (delta_trans_hat * cs_bearing + delta_strafe_hat * sn_bearing);
        sample->pose.v[1] +=
            (delta_trans_hat * sn_bearing - delta_strafe_hat * cs_bearing);
        sample->pose.v[2] += delta_rot_hat;
      }
    } break;
    case ODOM_MODEL_DIFF_CORRECTED: {
      // Implement sample_motion_odometry (Prob Rob p 136)
      double delta_rot1, delta_trans, delta_rot2;
      double delta_rot1_hat, delta_trans_hat, delta_rot2_hat;
      double delta_rot1_noise, delta_rot2_noise;

      // Avoid computing a bearing from two poses that are extremely near each
      // other (happens on in-place rotation).
      if (sqrt(ndata->delta.v[1] * ndata->delta.v[1] +
               ndata->delta.v[0] * ndata->delta.v[0]) < 0.01)
        delta_rot1 = 0.0;
      else
        delta_rot1 = angle_diff(atan2(ndata->delta.v[1], ndata->delta.v[0]),
                                old_pose.v[2]);
      delta_trans = sqrt(ndata->delta.v[0] * ndata->delta.v[0] +
                         ndata->delta.v[1] * ndata->delta.v[1]);
      delta_rot2 = angle_diff(ndata->delta.v[2], delta_rot1);

      // We want to treat backward and forward motion symmetrically for the
      // noise model to be applied below.  The standard model seems to assume
      // forward motion.
      delta_rot1_noise = std::min(fabs(angle_diff(delta_rot1, 0.0)),
                                  fabs(angle_diff(delta_rot1, M_PI)));
      delta_rot2_noise = std::min(fabs(angle_diff(delta_rot2, 0.0)),
                                  fabs(angle_diff(delta_rot2, M_PI)));

      for (int i = 0; i < set->samples.size(); i++) {
        pf::SamplePtr sample = set->samples[i];

        // Sample pose differences
        delta_rot1_hat = angle_diff(
            delta_rot1, pf::RandomGaussian(
                            sqrt(alpha1_ * delta_rot1_noise * delta_rot1_noise +
                                 alpha2_ * delta_trans * delta_trans)));
        delta_trans_hat =
            delta_trans - pf::RandomGaussian(sqrt(
                              alpha3_ * delta_trans * delta_trans +
                              alpha4_ * delta_rot1_noise * delta_rot1_noise +
                              alpha4_ * delta_rot2_noise * delta_rot2_noise));
        delta_rot2_hat = angle_diff(
            delta_rot2, pf::RandomGaussian(
                            sqrt(alpha1_ * delta_rot2_noise * delta_rot2_noise +
                                 alpha2_ * delta_trans * delta_trans)));

        // Apply sampled update to particle pose
        sample->pose.v[0] +=
            delta_trans_hat * cos(sample->pose.v[2] + delta_rot1_hat);
        sample->pose.v[1] +=
            delta_trans_hat * sin(sample->pose.v[2] + delta_rot1_hat);
        sample->pose.v[2] += delta_rot1_hat + delta_rot2_hat;
      }
    } break;
  }
  return true;
}
