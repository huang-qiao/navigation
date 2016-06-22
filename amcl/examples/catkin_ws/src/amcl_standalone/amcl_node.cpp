/*
 *  Copyright (c) 2008, Willow Garage, Inc.
 *  All rights reserved.
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/* Author: Brian Gerkey */

#include <algorithm>
#include <vector>
#include <map>
#include <cmath>

#include <boost/bind.hpp>
#include <boost/thread/mutex.hpp>

// Signal handling
#include <signal.h>

#include "map/map.h"
#include "pf/pf.h"
#include "sensors/amcl_odom.h"
#include "sensors/amcl_laser.h"
#include "amcl_imp.hpp"
/*
#include "asusbot/navigation/amcl/map/map.h"
#include "asusbot/navigation/amcl/pf/pf.h"
#include "asusbot/navigation/amcl/sensors/amcl_odom.h"
#include "asusbot/navigation/amcl/sensors/amcl_laser.h"
#include "asusbot/navigation/amcl/amcl_imp.hpp"
*/
#include "ros/assert.h"

// roscpp
#include "ros/ros.h"

// Messages that I need
#include "sensor_msgs/LaserScan.h"
#include "geometry_msgs/PoseWithCovarianceStamped.h"
#include "geometry_msgs/PoseArray.h"
#include "geometry_msgs/Pose.h"
#include "nav_msgs/GetMap.h"
#include "nav_msgs/SetMap.h"
#include "std_srvs/Empty.h"

// For transform support
#include "tf/transform_broadcaster.h"
#include "tf/transform_listener.h"
#include "tf/message_filter.h"
#include "tf/tf.h"
#include "tf_conversions/tf_eigen.h"
#include "message_filters/subscriber.h"

// Dynamic_reconfigure
#include "dynamic_reconfigure/server.h"

// Allows AMCL to run from bag file
#include <rosbag/bag.h>
#include <rosbag/view.h>
#include <boost/foreach.hpp>

#include <memory>

#include <Eigen/Geometry>

#define NEW_UNIFORM_SAMPLING 1

using namespace amcl;


static const std::string scan_topic_ = "scan";

class AmclNode : public AmclImpCallback
{
  public:
    AmclNode();
    ~AmclNode();

    inline void onUpdatePointCloud(pf_sample_set_t* set) override
    {
      geometry_msgs::PoseArray cloud_msg;
      cloud_msg.header.stamp = ros::Time::now();
      cloud_msg.header.frame_id = global_frame_id_;
      cloud_msg.poses.resize(set->sample_count);
      for(int i=0;i<set->sample_count;i++)
      {
        tf::poseTFToMsg(tf::Pose(tf::createQuaternionFromYaw(set->samples[i].pose.v[2]),
                                 tf::Vector3(set->samples[i].pose.v[0],
                                           set->samples[i].pose.v[1], 0)),
                        cloud_msg.poses[i]);
      }
      particlecloud_pub_.publish(cloud_msg);
    }

    inline void onUpdatePose(pf_vector_t pf_pose_mean, pf_matrix_t cov, double time_sec) override
    {
      geometry_msgs::PoseWithCovarianceStamped p;
      // Fill in the header
      p.header.frame_id = global_frame_id_;
      p.header.stamp.fromSec(time_sec);
      // Copy in the pose
      p.pose.pose.position.x = pf_pose_mean.v[0];
      p.pose.pose.position.y = pf_pose_mean.v[1];
      tf::quaternionTFToMsg(tf::createQuaternionFromYaw(pf_pose_mean.v[2]),
                            p.pose.pose.orientation);
      // Copy in the covariance, converting from 3-D to 6-D
      for(int i=0; i<2; i++)
      {
        for(int j=0; j<2; j++)
        {
          // Report the overall filter covariance, rather than the
          // covariance for the highest-weight cluster
          //p.covariance[6*i+j] = hyps[max_weight_hyp].pf_pose_cov.m[i][j];
          p.pose.covariance[6*i+j] = cov.m[i][j];
        }
      }
      // Report the overall filter covariance, rather than the
      // covariance for the highest-weight cluster
      //p.covariance[6*5+5] = hyps[max_weight_hyp].pf_pose_cov.m[2][2];
      p.pose.covariance[6*5+5] = cov.m[2][2];

      // debug only
      printf("cov:\n");
      for(int i=0; i<6; i++) {
        for(int j=0; j<6; j++) {
          printf("%6.3f ", p.pose.covariance[6*i+j]);
        }
        puts("");
      }

      pose_pub_.publish(p);
      last_published_pose = p;

      ROS_DEBUG("New pose: %6.3f %6.3f %6.3f",
               pf_pose_mean.v[0],
               pf_pose_mean.v[1],
               pf_pose_mean.v[2]);
    }

    inline void onUpdateTf(pf_vector_t pf_pose_mean, double time_sec) override
    {
      // subtracting base to odom from map to base and send map to odom instead
      tf::Stamped<tf::Pose> odom_to_map;
      try
      {
        tf::Transform tmp_tf(tf::createQuaternionFromYaw(pf_pose_mean.v[2]),
                             tf::Vector3(pf_pose_mean.v[0],
                                         pf_pose_mean.v[1],
                                         0.0));
        ros::Time stamp;
        stamp.fromSec(time_sec);
        tf::Stamped<tf::Pose> tmp_tf_stamped (tmp_tf.inverse(),
                                              stamp,
                                              base_frame_id_);
        this->tf_->transformPose(odom_frame_id_,
                                 tmp_tf_stamped,
                                 odom_to_map);
      }
      catch(tf::TransformException)
      {
        ROS_DEBUG("Failed to subtract base to odom transform");
        return;
      }

      tf::Transform latest_tf = tf::Transform(
            tf::Quaternion(
              odom_to_map.getRotation()),tf::Point(odom_to_map.getOrigin()
            )
      );

      if (tf_broadcast_ == true)
      {
        // We want to send a transform that is good up until a
        // tolerance time so that odom can be used
        ros::Time stamp;
        stamp.fromSec(time_sec);
        ros::Time transform_expiration = (stamp + transform_tolerance_);
        tf::StampedTransform tmp_tf_stamped(latest_tf.inverse(),
                                            transform_expiration,
                                            global_frame_id_, odom_frame_id_);
        this->tfb_->sendTransform(tmp_tf_stamped);
      }
    }

    inline void onQueryInitPose(AmclImp* amcl) override
    {
      // no need to use amcl ptr, updatePoseFromServer already take care of it.
      updatePoseFromServer();
    }

    inline void onSaveCurrentPose(pf_vector_t pf_pose, pf_matrix_t pf_cov) override
    {
      ros::Time now = ros::Time::now();
      if((save_pose_period.toSec() > 0.0) &&
         (now - save_pose_last_time) >= save_pose_period) {
        // no need for pf_pose and pf_cov, since savePoseToServer already take
        // care of it
        this->savePoseToServer();
        save_pose_last_time = now;
      }
    }

    /**
     * @brief Uses TF and LaserScan messages from bag file to drive AMCL instead
     */
    void runFromBag(const std::string &in_bag_fn);

    int process();
    void savePoseToServer();

  private:
    tf::TransformBroadcaster* tfb_;

    // Use a child class to get access to tf2::Buffer class inside of tf_
    struct TransformListenerWrapper : public tf::TransformListener
    {
      inline tf2_ros::Buffer &getBuffer() {return tf2_buffer_;}
    };

    TransformListenerWrapper* tf_;


    tf::Transform latest_tf_;
    // Callbacks
    bool globalLocalizationCallback(std_srvs::Empty::Request& req,
                                    std_srvs::Empty::Response& res);
    bool nomotionUpdateCallback(std_srvs::Empty::Request& req,
                                    std_srvs::Empty::Response& res);
    bool setMapCallback(nav_msgs::SetMap::Request& req,
                        nav_msgs::SetMap::Response& res);

    void laserReceived(const sensor_msgs::LaserScanConstPtr& laser_scan);
    void initialPoseReceived(const geometry_msgs::PoseWithCovarianceStampedConstPtr& msg);
    void handleInitialPoseMessage(const geometry_msgs::PoseWithCovarianceStamped& msg);
    void mapReceived(const nav_msgs::OccupancyGridConstPtr& msg);

    void handleMapMessage(const nav_msgs::OccupancyGrid& msg);
    map_t* convertMap( const nav_msgs::OccupancyGrid& map_msg );
    void updatePoseFromServer();

    double getYaw(tf::Pose& t);

    //parameter for what odom to use
    std::string odom_frame_id_;

    //paramater to store latest odom pose
    //tf::Stamped<tf::Pose> latest_odom_pose_;
    //bool latest_tf_valid_;
    bool tf_broadcast_;

    //parameter for what base to use
    std::string base_frame_id_;
    std::string global_frame_id_;

    bool use_map_topic_;

    ros::Duration gui_publish_period;
    ros::Time save_pose_last_time;
    ros::Duration save_pose_period;

    geometry_msgs::PoseWithCovarianceStamped last_published_pose;

    message_filters::Subscriber<sensor_msgs::LaserScan>* laser_scan_sub_;
    tf::MessageFilter<sensor_msgs::LaserScan>* laser_scan_filter_;
    ros::Subscriber initial_pose_sub_;

    ros::Duration cloud_pub_interval;
    ros::Time last_cloud_pub_time;

    // For slowing play-back when reading directly from a bag file
    ros::WallDuration bag_scan_period_;

    void requestMap();

    // Helper to get odometric pose from transform system
    bool getOdomPose(tf::Stamped<tf::Pose>& pose,
                     double& x, double& y, double& yaw,
                     const ros::Time& t, const std::string& f);

    //time for tolerance on the published transform,
    //basically defines how long a map->odom transform is good for
    ros::Duration transform_tolerance_;

    ros::NodeHandle nh_;
    ros::NodeHandle private_nh_;
    ros::Publisher pose_pub_;
    ros::Publisher particlecloud_pub_;
    ros::ServiceServer global_loc_srv_;
    ros::ServiceServer nomotion_update_srv_; //to let amcl update samples without requiring motion
    ros::ServiceServer set_map_srv_;
    ros::Subscriber initial_pose_sub_old_;
    ros::Subscriber map_sub_;

    ros::Timer check_laser_timer_;

    ros::Time last_laser_received_ts_;
    ros::Duration laser_check_interval_;
    void checkLaserReceived(const ros::TimerEvent& event);

    std::shared_ptr<AmclImp> amcl_;
};


#define USAGE "USAGE: amcl"

boost::shared_ptr<AmclNode> amcl_node_ptr;

void sigintHandler(int sig)
{
  // Save latest pose as we're shutting down.
  amcl_node_ptr->savePoseToServer();
  ros::shutdown();
}

int
main(int argc, char** argv)
{
  ros::init(argc, argv, "amcl");
  ros::NodeHandle nh;

  // Override default sigint handler
  signal(SIGINT, sigintHandler);

  // Make our node available to sigintHandler
  amcl_node_ptr.reset(new AmclNode());

  if (argc == 1)
  {
    // run using ROS input
    ros::spin();
  }
  else if ((argc == 3) && (std::string(argv[1]) == "--run-from-bag"))
  {
    amcl_node_ptr->runFromBag(argv[2]);
  }

  // Without this, our boost locks are not shut down nicely
  amcl_node_ptr.reset();

  // To quote Morgan, Hooray!
  return(0);
}

AmclNode::AmclNode() :
  tf_broadcast_(true),
  private_nh_("~")
{
  amcl_ = std::make_shared<AmclImp>();
  amcl_->registerCallback(this);

  // Grab params off the param server
  private_nh_.param("use_map_topic", use_map_topic_, false);

  double tmp;
  private_nh_.param("gui_publish_rate", tmp, -1.0);
  gui_publish_period = ros::Duration(1.0/tmp);
  private_nh_.param("save_pose_rate", tmp, 0.5);
  save_pose_period = ros::Duration(1.0/tmp);

  private_nh_.param("odom_frame_id", odom_frame_id_, std::string("odom"));
  private_nh_.param("base_frame_id", base_frame_id_, std::string("base_link"));
  private_nh_.param("global_frame_id", global_frame_id_, std::string("map"));
  double tmp_tol;
  private_nh_.param("transform_tolerance", tmp_tol, 0.1);

  transform_tolerance_.fromSec(tmp_tol);

  {
    double bag_scan_period;
    private_nh_.param("bag_scan_period", bag_scan_period, -1.0);
    bag_scan_period_.fromSec(bag_scan_period);
  }

  updatePoseFromServer();

  cloud_pub_interval.fromSec(1.0);
  tfb_ = new tf::TransformBroadcaster();
  tf_ = new TransformListenerWrapper();

  pose_pub_ = nh_.advertise<geometry_msgs::PoseWithCovarianceStamped>("amcl_pose", 2, true);
  particlecloud_pub_ = nh_.advertise<geometry_msgs::PoseArray>("particlecloud", 2, true);
  global_loc_srv_ = nh_.advertiseService("global_localization",
                                         &AmclNode::globalLocalizationCallback,
                                         this);
  nomotion_update_srv_= nh_.advertiseService("request_nomotion_update", &AmclNode::nomotionUpdateCallback, this);
  set_map_srv_= nh_.advertiseService("set_map", &AmclNode::setMapCallback, this);

  laser_scan_sub_ = new message_filters::Subscriber<sensor_msgs::LaserScan>(nh_, scan_topic_, 100);
  laser_scan_filter_ =
          new tf::MessageFilter<sensor_msgs::LaserScan>(*laser_scan_sub_,
                                                        *tf_,
                                                        odom_frame_id_,
                                                        100);
  laser_scan_filter_->registerCallback(boost::bind(&AmclNode::laserReceived,
                                                   this, _1));
  initial_pose_sub_ = nh_.subscribe("initialpose", 2, &AmclNode::initialPoseReceived, this);

  if(use_map_topic_) {
    map_sub_ = nh_.subscribe("map", 1, &AmclNode::mapReceived, this);
    ROS_INFO("Subscribed to map topic.");
  } else {
    requestMap();
  }

  // 15s timer to warn on lack of receipt of laser scans, #5209
  laser_check_interval_ = ros::Duration(15.0);
  check_laser_timer_ = nh_.createTimer(laser_check_interval_,
                                       boost::bind(&AmclNode::checkLaserReceived, this, _1));
}

void AmclNode::runFromBag(const std::string &in_bag_fn)
{
  rosbag::Bag bag;
  bag.open(in_bag_fn, rosbag::bagmode::Read);
  std::vector<std::string> topics;
  topics.push_back(std::string("tf"));
  std::string scan_topic_name = "base_scan"; // TODO determine what topic this actually is from ROS
  topics.push_back(scan_topic_name);
  rosbag::View view(bag, rosbag::TopicQuery(topics));

  ros::Publisher laser_pub = nh_.advertise<sensor_msgs::LaserScan>(scan_topic_name, 100);
  ros::Publisher tf_pub = nh_.advertise<tf2_msgs::TFMessage>("/tf", 100);

  // Sleep for a second to let all subscribers connect
  ros::WallDuration(1.0).sleep();

  ros::WallTime start(ros::WallTime::now());

  // Wait for map
  while (ros::ok())
  {
    {
      if (amcl_->getMap())
      {
        ROS_INFO("Map is ready");
        break;
      }
    }
    ROS_INFO("Waiting for map...");
    ros::getGlobalCallbackQueue()->callAvailable(ros::WallDuration(1.0));
  }

  BOOST_FOREACH(rosbag::MessageInstance const msg, view)
  {
    if (!ros::ok())
    {
      break;
    }

    // Process any ros messages or callbacks at this point
    ros::getGlobalCallbackQueue()->callAvailable(ros::WallDuration());

    tf2_msgs::TFMessage::ConstPtr tf_msg = msg.instantiate<tf2_msgs::TFMessage>();
    if (tf_msg != NULL)
    {
      tf_pub.publish(msg);
      for (size_t ii=0; ii<tf_msg->transforms.size(); ++ii)
      {
        tf_->getBuffer().setTransform(tf_msg->transforms[ii], "rosbag_authority");
      }
      continue;
    }

    sensor_msgs::LaserScan::ConstPtr base_scan = msg.instantiate<sensor_msgs::LaserScan>();
    if (base_scan != NULL)
    {
      laser_pub.publish(msg);
      laser_scan_filter_->add(base_scan);
      if (bag_scan_period_ > ros::WallDuration(0))
      {
        bag_scan_period_.sleep();
      }
      continue;
    }

    ROS_WARN_STREAM("Unsupported message type" << msg.getTopic());
  }

  bag.close();

  double runtime = (ros::WallTime::now() - start).toSec();
  ROS_INFO("Bag complete, took %.1f seconds to process, shutting down", runtime);

  const geometry_msgs::Quaternion & q(last_published_pose.pose.pose.orientation);
  double yaw, pitch, roll;
  tf::Matrix3x3(tf::Quaternion(q.x, q.y, q.z, q.w)).getEulerYPR(yaw,pitch,roll);
  ROS_INFO("Final location %.3f, %.3f, %.3f with stamp=%f",
            last_published_pose.pose.pose.position.x,
            last_published_pose.pose.pose.position.y,
            yaw, last_published_pose.header.stamp.toSec()
            );

  ros::shutdown();
}


void AmclNode::savePoseToServer()
{
  // We need to apply the last transform to the latest odom pose to get
  // the latest map pose to store.  We'll take the covariance from
  // last_published_pose.
  pf_vector_t pf_odom_pose = amcl_->getLatestOdomPose();
  tf::Stamped<tf::Pose> latest_odom_pose;
  latest_odom_pose.setOrigin(
        tf::Vector3(pf_odom_pose.v[0], pf_odom_pose.v[1], 0));
  tf::Matrix3x3 matRot;
  matRot.setRPY(0,0,pf_odom_pose.v[2]);
  latest_odom_pose.setBasis(matRot);
  tf::Pose map_pose = latest_tf_.inverse() * latest_odom_pose;
  double yaw,pitch,roll;
  map_pose.getBasis().getEulerYPR(yaw, pitch, roll);

  ROS_DEBUG("Saving pose to server. x: %.3f, y: %.3f", map_pose.getOrigin().x(), map_pose.getOrigin().y() );

  private_nh_.setParam("initial_pose_x", map_pose.getOrigin().x());
  private_nh_.setParam("initial_pose_y", map_pose.getOrigin().y());
  private_nh_.setParam("initial_pose_a", yaw);
  private_nh_.setParam("initial_cov_xx",
                                  last_published_pose.pose.covariance[6*0+0]);
  private_nh_.setParam("initial_cov_yy",
                                  last_published_pose.pose.covariance[6*1+1]);
  private_nh_.setParam("initial_cov_aa",
                                  last_published_pose.pose.covariance[6*5+5]);

}

void AmclNode::updatePoseFromServer()
{
  double init_pose[3];
  double init_cov[3];
  init_pose[0] = 0.0;
  init_pose[1] = 0.0;
  init_pose[2] = 0.0;
  init_cov[0] = 0.5 * 0.5;
  init_cov[1] = 0.5 * 0.5;
  init_cov[2] = (M_PI/12.0) * (M_PI/12.0);
  amcl_->setPose(init_pose, init_cov);
  // Check for NAN on input from param server, #5239
  double tmp_pos;
  private_nh_.param("initial_pose_x", tmp_pos, init_pose[0]);
  if(!std::isnan(tmp_pos))
    init_pose[0] = tmp_pos;
  else
    ROS_WARN("ignoring NAN in initial pose X position");
  private_nh_.param("initial_pose_y", tmp_pos, init_pose[1]);
  if(!std::isnan(tmp_pos))
    init_pose[1] = tmp_pos;
  else
    ROS_WARN("ignoring NAN in initial pose Y position");
  private_nh_.param("initial_pose_a", tmp_pos, init_pose[2]);
  if(!std::isnan(tmp_pos))
    init_pose[2] = tmp_pos;
  else
    ROS_WARN("ignoring NAN in initial pose Yaw");
  private_nh_.param("initial_cov_xx", tmp_pos, init_cov[0]);
  if(!std::isnan(tmp_pos))
    init_cov[0] =tmp_pos;
  else
    ROS_WARN("ignoring NAN in initial covariance XX");
  private_nh_.param("initial_cov_yy", tmp_pos, init_cov[1]);
  if(!std::isnan(tmp_pos))
    init_cov[1] = tmp_pos;
  else
    ROS_WARN("ignoring NAN in initial covariance YY");
  private_nh_.param("initial_cov_aa", tmp_pos, init_cov[2]);
  if(!std::isnan(tmp_pos))
    init_cov[2] = tmp_pos;
  else
    ROS_WARN("ignoring NAN in initial covariance AA");

  amcl_->setPose(init_pose, init_cov);

}

void
AmclNode::checkLaserReceived(const ros::TimerEvent& event)
{

  ros::Duration d = ros::Time::now() - last_laser_received_ts_;
  if(d > laser_check_interval_)
  {
    ROS_WARN("No laser scan received (and thus no pose updates have been published) for %f seconds.  Verify that data is being published on the %s topic.",
             d.toSec(),
             ros::names::resolve(scan_topic_).c_str());
  }

}

void
AmclNode::requestMap()
{

  // get map via RPC
  nav_msgs::GetMap::Request  req;
  nav_msgs::GetMap::Response resp;
  ROS_INFO("Requesting the map...");
  while(!ros::service::call("static_map", req, resp))
  {
    ROS_WARN("Request for map failed; trying again...");
    ros::Duration d(0.5);
    d.sleep();
  }

  handleMapMessage( resp.map );

}

void
AmclNode::mapReceived(const nav_msgs::OccupancyGridConstPtr& msg)
{
  handleMapMessage( *msg );
}

void
AmclNode::handleMapMessage(const nav_msgs::OccupancyGrid& msg)
{
  ROS_INFO("Received a %d X %d map @ %.3f m/pix\n",
           msg.info.width,
           msg.info.height,
           msg.info.resolution);

  amcl_->freeMapDependentMemory();
  // Clear queued laser objects because they hold pointers to the existing map, #5202.
  amcl_->initMap(convertMap(msg));
}

/**
 * Convert an OccupancyGrid map message into the internal
 * representation.  This allocates a map_t and returns it.
 */
map_t*
AmclNode::convertMap( const nav_msgs::OccupancyGrid& map_msg )
{

  map_t* map = map_alloc();
  ROS_ASSERT(map);

  map->size_x = map_msg.info.width;
  map->size_y = map_msg.info.height;
  map->scale = map_msg.info.resolution;
  map->origin_x = map_msg.info.origin.position.x + (map->size_x / 2) * map->scale;
  map->origin_y = map_msg.info.origin.position.y + (map->size_y / 2) * map->scale;
  // Convert to player format
  map->cells = (map_cell_t*)malloc(sizeof(map_cell_t)*map->size_x*map->size_y);
  ROS_ASSERT(map->cells);
  for(int i=0;i<map->size_x * map->size_y;i++)
  {
    if(map_msg.data[i] == 0)
      map->cells[i].occ_state = -1;
    else if(map_msg.data[i] == 100)
      map->cells[i].occ_state = +1;
    else
      map->cells[i].occ_state = 0;
  }

  return map;
}

AmclNode::~AmclNode()
{
  amcl_->freeMapDependentMemory();
  delete laser_scan_filter_;
  delete laser_scan_sub_;
  delete tfb_;
  delete tf_;
  // TODO: delete everything allocated in constructor
}

bool
AmclNode::getOdomPose(tf::Stamped<tf::Pose>& odom_pose,
                      double& x, double& y, double& yaw,
                      const ros::Time& t, const std::string& f)
{

  // Get the robot's pose
  tf::Stamped<tf::Pose> ident (tf::Transform(tf::createIdentityQuaternion(),
                                           tf::Vector3(0,0,0)), t, f);
  try
  {
    this->tf_->transformPose(odom_frame_id_, ident, odom_pose);
  }
  catch(tf::TransformException e)
  {
    ROS_WARN("Failed to compute odom pose, skipping scan (%s)", e.what());
    return false;
  }
  x = odom_pose.getOrigin().x();
  y = odom_pose.getOrigin().y();
  double pitch,roll;
  odom_pose.getBasis().getEulerYPR(yaw, pitch, roll);

  return true;
}



bool
AmclNode::globalLocalizationCallback(std_srvs::Empty::Request& req,
                                     std_srvs::Empty::Response& res)
{
  if( amcl_->getMap() == NULL ) {
    return true;
  }
  return amcl_->doGlobalLocalization();
}

// force nomotion updates (amcl updating without requiring motion)
bool
AmclNode::nomotionUpdateCallback(std_srvs::Empty::Request& req,
                                     std_srvs::Empty::Response& res)
{
  return amcl_->doNoMotionUpdate();
}

bool
AmclNode::setMapCallback(nav_msgs::SetMap::Request& req,
                         nav_msgs::SetMap::Response& res)
{
  handleMapMessage(req.map);
  handleInitialPoseMessage(req.initial_pose);
  res.success = true;
  return true;
}

void
AmclNode::laserReceived(const sensor_msgs::LaserScanConstPtr& laser_scan)
{
  last_laser_received_ts_ = ros::Time::now();
  if( amcl_->getMap() == NULL ) {
    return;
  }
  tf::Stamped<tf::Pose> ident (tf::Transform(tf::createIdentityQuaternion(), tf::Vector3(0,0,0)), ros::Time(), laser_scan->header.frame_id);
  tf::Stamped<tf::Pose> laser_pose;
  try
  {
    this->tf_->transformPose(base_frame_id_, ident, laser_pose);
  }
  catch(tf::TransformException& e)
  {
    ROS_ERROR("Couldn't transform from %s to %s, "
              "even though the message notifier is in use", laser_scan->header.frame_id.c_str(), base_frame_id_.c_str());
    return;
  }

  Eigen::Isometry3d se3_laser_pose;
  tf::poseTFToEigen(laser_pose, se3_laser_pose);

  // Where was the robot when this scan was taken?
  pf_vector_t pose;
  tf::Stamped<tf::Pose> latest_odom_pose;
  if(!getOdomPose(latest_odom_pose, pose.v[0], pose.v[1], pose.v[2], laser_scan->header.stamp, base_frame_id_))
  {
    ROS_ERROR("Couldn't determine robot's pose associated with laser scan");
    return;
  }

  amcl_->processSensorPose(laser_scan->header.frame_id, se3_laser_pose, pose);

  // If the robot has moved, update the filter
  if(amcl_->isLaserUpdate(laser_scan->header.frame_id))
  {
    AMCLLaserData ldata;
    //ldata.sensor = lasers_[laser_index];
    ldata.sensor = amcl_->getLaser(laser_scan->header.frame_id);
    ldata.range_count = laser_scan->ranges.size();
    ldata.time = laser_scan->header.stamp.toSec();

    // To account for lasers that are mounted upside-down, we determine the
    // min, max, and increment angles of the laser in the base frame.
    //
    // Construct min and max angles of laser, in the base_link frame.
    tf::Quaternion q;
    q.setRPY(0.0, 0.0, laser_scan->angle_min);
    tf::Stamped<tf::Quaternion> min_q(q, laser_scan->header.stamp,
                                      laser_scan->header.frame_id);
    q.setRPY(0.0, 0.0, laser_scan->angle_min + laser_scan->angle_increment);
    tf::Stamped<tf::Quaternion> inc_q(q, laser_scan->header.stamp,
                                      laser_scan->header.frame_id);
    try
    {
      tf_->transformQuaternion(base_frame_id_, min_q, min_q);
      tf_->transformQuaternion(base_frame_id_, inc_q, inc_q);
    }
    catch(tf::TransformException& e)
    {
      ROS_WARN("Unable to transform min/max laser angles into base frame: %s",
               e.what());
      return;
    }

    double angle_min = tf::getYaw(min_q);
    double angle_increment = tf::getYaw(inc_q) - angle_min;

    // wrapping angle to [-pi .. pi]
    angle_increment = fmod(angle_increment + 5*M_PI, 2*M_PI) - M_PI;

    ROS_DEBUG("Laser angles in base frame: min: %.3f inc: %.3f", angle_min, angle_increment);

    // Apply range min/max thresholds, if the user supplied them
    if (amcl_->getLaserMaxRange() > 0.0)
      ldata.range_max = std::min(laser_scan->range_max, (float)amcl_->getLaserMaxRange());
    else
      ldata.range_max = laser_scan->range_max;
    double range_min;
    if(amcl_->getLaserMinRange() > 0.0)
      range_min = std::max(laser_scan->range_min, (float)amcl_->getLaserMinRange());
    else
      range_min = laser_scan->range_min;
    // The AMCLLaserData destructor will free this memory
    ldata.ranges = new double[ldata.range_count][2];
    ROS_ASSERT(ldata.ranges);
    for(int i=0;i<ldata.range_count;i++)
    {
      // amcl doesn't (yet) have a concept of min range.  So we'll map short
      // readings to max range.
      if(laser_scan->ranges[i] <= range_min)
        ldata.ranges[i][0] = ldata.range_max;
      else
        ldata.ranges[i][0] = laser_scan->ranges[i];
      // Compute bearing
      ldata.ranges[i][1] = angle_min + (i * angle_increment);
    }

    amcl_->processLaserData(laser_scan->header.frame_id, &ldata, pose);
  }

}

double
AmclNode::getYaw(tf::Pose& t)
{
  double yaw, pitch, roll;
  t.getBasis().getEulerYPR(yaw,pitch,roll);
  return yaw;
}

void
AmclNode::initialPoseReceived(const geometry_msgs::PoseWithCovarianceStampedConstPtr& msg)
{
  handleInitialPoseMessage(*msg);
}

void
AmclNode::handleInitialPoseMessage(const geometry_msgs::PoseWithCovarianceStamped& msg)
{
  if(msg.header.frame_id == "")
  {
    // This should be removed at some point
    ROS_WARN("Received initial pose with empty frame_id.  You should always supply a frame_id.");
  }
  // We only accept initial pose estimates in the global frame, #5148.
  else if(tf_->resolve(msg.header.frame_id) != tf_->resolve(global_frame_id_))
  {
    ROS_WARN("Ignoring initial pose in frame \"%s\"; initial poses must be in the global frame, \"%s\"",
             msg.header.frame_id.c_str(),
             global_frame_id_.c_str());
    return;
  }

  // In case the client sent us a pose estimate in the past, integrate the
  // intervening odometric change.
  tf::StampedTransform tx_odom;
  try
  {
    ros::Time now = ros::Time::now();
    // wait a little for the latest tf to become available
    tf_->waitForTransform(base_frame_id_, msg.header.stamp,
                         base_frame_id_, now,
                         odom_frame_id_, ros::Duration(0.5));
    tf_->lookupTransform(base_frame_id_, msg.header.stamp,
                         base_frame_id_, now,
                         odom_frame_id_, tx_odom);
  }
  catch(tf::TransformException e)
  {
    // If we've never sent a transform, then this is normal, because the
    // global_frame_id_ frame doesn't exist.  We only care about in-time
    // transformation for on-the-move pose-setting, so ignoring this
    // startup condition doesn't really cost us anything.
    if(amcl_->sentFirstTransform())
      ROS_WARN("Failed to transform initial pose in time (%s)", e.what());
    tx_odom.setIdentity();
  }

  tf::Pose pose_old, pose_new;
  tf::poseMsgToTF(msg.pose.pose, pose_old);
  pose_new = pose_old * tx_odom;

  // Transform into the global frame

  ROS_INFO("Setting pose (%.6f): %.3f %.3f %.3f",
           ros::Time::now().toSec(),
           pose_new.getOrigin().x(),
           pose_new.getOrigin().y(),
           getYaw(pose_new));
  // Re-initialize the filter
  pf_vector_t pf_init_pose_mean = pf_vector_zero();
  pf_init_pose_mean.v[0] = pose_new.getOrigin().x();
  pf_init_pose_mean.v[1] = pose_new.getOrigin().y();
  pf_init_pose_mean.v[2] = getYaw(pose_new);
  pf_matrix_t pf_init_pose_cov = pf_matrix_zero();
  // Copy in the covariance, converting from 6-D to 3-D
  for(int i=0; i<2; i++)
  {
    for(int j=0; j<2; j++)
    {
      pf_init_pose_cov.m[i][j] = msg.pose.covariance[6*i+j];
    }
  }
  pf_init_pose_cov.m[2][2] = msg.pose.covariance[6*5+5];

  amcl_->setInitPose(pf_init_pose_mean, pf_init_pose_cov);
}

