#pragma once

#include "amcl_sensor.hpp"
#include "../pf/pf_pdf.hpp"

namespace amcl {

struct AMCLOdomData;
using AMCLOdomDataPtr = std::shared_ptr<AMCLOdomData>;

class AMCLOdom;
using AMCLOdomPtr = std::shared_ptr<AMCLOdom>;

typedef enum
{
  ODOM_MODEL_DIFF,
  ODOM_MODEL_OMNI,
  ODOM_MODEL_DIFF_CORRECTED,
  ODOM_MODEL_OMNI_CORRECTED
} OdomModelType;

// Odometric sensor data
struct AMCLOdomData : public AMCLSensorData
{
  // Odometric pose
  pf::Pose pose;

  // Change in odometric pose
  pf::Pose delta;
};


// Odometric sensor model
class AMCLOdom : public AMCLSensor
{
public:
  // Default constructor
  AMCLOdom();

  void setModelDiff(double alpha1,
                    double alpha2,
                    double alpha3,
                    double alpha4);

  void setModelOmni(double alpha1,
                    double alpha2,
                    double alpha3,
                    double alpha4,
                    double alpha5);

  void setModel(OdomModelType type,
                double alpha1,
                double alpha2,
                double alpha3,
                double alpha4,
                double alpha5 = 0 );

  // Update the filter based on the action model.  Returns true if the filter
  // has been updated.
  virtual bool updateAction(pf::ParticleFilterPtr pf, AMCLSensorDataPtr data);

private:
  // Current data timestamp
  double time_;
  
  // Model type
  OdomModelType model_type_;

  // Drift parameters
  double alpha1_, alpha2_, alpha3_, alpha4_, alpha5_;
};


} // namespace amcl
