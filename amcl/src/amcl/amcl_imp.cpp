#include "amcl/amcl_imp.hpp"

#include <vector>
using namespace amcl;
std::vector<std::pair<int, int> > AmclImp::free_space_indices;

pf_vector_t AmclImp::uniformPoseGenerator(void* arg)
{
   map_t* map = (map_t*)arg;
#if NEW_UNIFORM_SAMPLING
   unsigned int rand_index = drand48() * free_space_indices.size();
   std::pair<int,int> free_point = free_space_indices[rand_index];
   pf_vector_t p;
   p.v[0] = MAP_WXGX(map, free_point.first);
   p.v[1] = MAP_WYGY(map, free_point.second);
   p.v[2] = drand48() * 2 * M_PI - M_PI;
#else
   double min_x, max_x, min_y, max_y;

   min_x = (map->size_x * map->scale)/2.0 - map->origin_x;
   max_x = (map->size_x * map->scale)/2.0 + map->origin_x;
   min_y = (map->size_y * map->scale)/2.0 - map->origin_y;
   max_y = (map->size_y * map->scale)/2.0 + map->origin_y;

   pf_vector_t p;

   DLOG(INFO) << "Generating new uniform sample";
   for(;;)
   {
     p.v[0] = min_x + drand48() * (max_x - min_x);
     p.v[1] = min_y + drand48() * (max_y - min_y);
     p.v[2] = drand48() * 2 * M_PI - M_PI;
     // Check that it's a free cell
     int i,j;
     i = MAP_GXWX(map, p.v[0]);
     j = MAP_GYWY(map, p.v[1]);
     if(MAP_VALID(map,i,j) && (map->cells[MAP_INDEX(map,i,j)].occ_state == -1))
       break;
     }
#endif
  return p;
}

AmclImp::AmclImp()
    : sent_first_transform_(false),
      map_(NULL),
      pf_(NULL),
      resample_count_(0),
      odom_(NULL),
      laser_(NULL),
      initial_pose_hyp_(NULL),
      first_map_received_(false),
      first_reconfigure_call_(true) {
  // Grab params off the param server
  first_map_only_ = false;

  laser_min_range_ = -1.0;
  laser_max_range_ = -1.0;
  max_beams_ = 30;
  min_particles_ = 100;
  max_particles_ = 5000;
  pf_err_ = 0.01;
  pf_z_ = 0.99;
  alpha1_ = 0.2;
  alpha2_ = 0.2;
  alpha3_ = 0.2;
  alpha4_ = 0.2;
  alpha5_ = 0.2;

  do_beamskip_ = false;
  beam_skip_distance_ = 0.5;
  beam_skip_threshold_ = 0.3;
  beam_skip_error_threshold_ = 0.9;

  z_hit_ = 0.95;
  z_short_ = 0.1;
  z_max_ = 0.05;
  z_rand_ = 0.05;
  sigma_hit_ = 0.2;
  lambda_short_ = 0.1;
  laser_likelihood_max_dist_ = 2.0;
  // laser_model_type_ = LASER_MODEL_BEAM;
  laser_model_type_ = LASER_MODEL_LIKELIHOOD_FIELD;
  // laser_model_type_ = LASER_MODEL_LIKELIHOOD_FIELD_PROB;
  odom_model_type_ = ODOM_MODEL_DIFF;
  // odom_model_type_ = ODOM_MODEL_OMNI;
  // odom_model_type_ = ODOM_MODEL_DIFF_CORRECTED;
  // odom_model_type_ = ODOM_MODEL_OMNI_CORRECTED;

  d_thresh_ = 0.2;
  a_thresh_ = M_PI / 6.0;
  resample_interval_ = 2;
  alpha_slow_ = 0.001;
  alpha_fast_ = 0.1;

  m_force_update = false;
}

void AmclImp::processSensorPose(std::string frame_id, Eigen::Isometry3d laser_pose, pf_vector_t odom_pose)
{
  int laser_index = -1;
  // Do we have the base->base_laser Tx yet?
  if(frame_to_laser_.find(frame_id) == frame_to_laser_.end())
  {
    printf("Setting up laser %d (frame_id=%s)\n", (int)frame_to_laser_.size(), frame_id.c_str());
    lasers_.push_back(new AMCLLaser(*laser_));
    lasers_update_.push_back(true);
    laser_index = frame_to_laser_.size();
    pf_vector_t laser_pose_v;
    laser_pose_v.v[0] = laser_pose.translation().matrix()(0);//getOrigin().x();
    laser_pose_v.v[1] = laser_pose.translation().matrix()(1); //getOrigin().y();
    // laser mounting angle gets computed later -> set to 0 here!
    laser_pose_v.v[2] = 0;
    lasers_[laser_index]->SetLaserPose(laser_pose_v);
    printf("Received laser's pose wrt robot: %.3f %.3f %.3f\n",
              laser_pose_v.v[0],
              laser_pose_v.v[1],
              laser_pose_v.v[2]);

    frame_to_laser_[frame_id] = laser_index;
  } else {
    // we have the laser pose, retrieve laser index
    laser_index = frame_to_laser_[frame_id];
  }

  pf_vector_t delta = pf_vector_zero();

  if(pf_init_)
  {

    // Compute change in pose
    //delta = pf_vector_coord_sub(pose, pf_odom_pose_);
    delta.v[0] = odom_pose.v[0] - pf_odom_pose_.v[0];
    delta.v[1] = odom_pose.v[1] - pf_odom_pose_.v[1];
    delta.v[2] = angle_diff(odom_pose.v[2], pf_odom_pose_.v[2]);

    // See if we should update the filter
    bool update = fabs(delta.v[0]) > d_thresh_ ||
                  fabs(delta.v[1]) > d_thresh_ ||
                  fabs(delta.v[2]) > a_thresh_;
    update = update || m_force_update;
    m_force_update=false;

    // Set the laser update flags
    if(update)
      for(unsigned int i=0; i < lasers_update_.size(); i++)
        lasers_update_[i] = true;
  }

  force_publication_ = false;
  if(!pf_init_)
  {

    // Pose at last filter update
    pf_odom_pose_ = odom_pose;

    // Filter is now initialized
    pf_init_ = true;

    // Should update sensor data
    for(unsigned int i=0; i < lasers_update_.size(); i++)
      lasers_update_[i] = true;

    force_publication_ = true;

    resample_count_ = 0;
  }
  // If the robot has moved, update the filter
  else if(pf_init_ && lasers_update_[laser_index])
  {

    printf("pose\n");
    pf_vector_fprintf(odom_pose, stdout, "%.3f");

    AMCLOdomData odata;
    odata.pose = odom_pose;
    // HACK
    // Modify the delta in the action data so the filter gets
    // updated correctly
    odata.delta = delta;

    // Use the action data to update the filter
    odom_->UpdateAction(pf_, (AMCLSensorData*)&odata);

    // Pose at last filter update
    //this->pf_odom_pose = pose;
    pf_sample_set_t* set = pf_->sets + pf_->current_set;
    onUpdatePointCloud(set);

  }
}

void AmclImp::processLaserData(std::string frame_id, AMCLLaserData* ldata, pf_vector_t odom_pose)
{
  bool resampled = false;
  lasers_[frame_to_laser_[frame_id]]->UpdateSensor(pf_, (AMCLSensorData*)ldata);
  lasers_update_[frame_to_laser_[frame_id]] = false;

  pf_odom_pose_ = odom_pose;

  // Resample the particles
  if(!(++resample_count_ % resample_interval_))
  {
    pf_update_resample(pf_);
    resampled = true;
  }

  pf_sample_set_t* set = pf_->sets + pf_->current_set;
  printf("Num samples: %d\n", set->sample_count);

  // Publish the resulting cloud
  // TODO: set maximum rate for publishing
  if (!m_force_update) {
    onUpdatePointCloud(set);
  }

  if(resampled || force_publication_)
  {
    // Read out the current hypotheses
    double max_weight = 0.0;
    int max_weight_hyp = -1;
    std::vector<amcl_hyp_t> hyps;
    hyps.resize(pf_->sets[pf_->current_set].cluster_count);
    for(int hyp_count = 0;
        hyp_count < pf_->sets[pf_->current_set].cluster_count; hyp_count++)
    {
      double weight;
      pf_vector_t pose_mean;
      pf_matrix_t pose_cov;
      if (!pf_get_cluster_stats(pf_, hyp_count, &weight, &pose_mean, &pose_cov))
      {
        printf("Couldn't get stats on cluster %d", hyp_count);
        break;
      }

      hyps[hyp_count].weight = weight;
      hyps[hyp_count].pf_pose_mean = pose_mean;
      hyps[hyp_count].pf_pose_cov = pose_cov;

      if(hyps[hyp_count].weight > max_weight)
      {
        max_weight = hyps[hyp_count].weight;
        max_weight_hyp = hyp_count;
      }
    }

    if(max_weight > 0.0)
    {
      printf("Max weight pose: %.3f %.3f %.3f\n",
                hyps[max_weight_hyp].pf_pose_mean.v[0],
                hyps[max_weight_hyp].pf_pose_mean.v[1],
                hyps[max_weight_hyp].pf_pose_mean.v[2]);

      puts("");
      pf_matrix_fprintf(hyps[max_weight_hyp].pf_pose_cov, stdout, "%6.3f");
      puts("");

      pf_sample_set_t* set = pf_->sets + pf_->current_set;

      onUpdatePose(hyps[max_weight_hyp].pf_pose_mean, set->cov, ldata->time);

      onUpdateTf(hyps[max_weight_hyp].pf_pose_mean, ldata->time);

      latest_pose_ = hyps[max_weight_hyp].pf_pose_mean;
      latest_cov_ = set->cov;
      latest_pose_valid_ = true;
    }
    else
    {
      printf("No pose!");
    }
  } else if(latest_pose_valid_) {
    // Nothing changed, so we'll just republish the last transform, to keep
    // everybody happy.
    onUpdateTf(latest_pose_, ldata->time);

    // Is it time to save our last pose to the param server
    onSavePose();
  }
}
