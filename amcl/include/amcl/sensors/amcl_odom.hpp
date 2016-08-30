#pragma once

#include "amcl_sensor.hpp"
#include "../pf/pf_pdf.hpp"

namespace amcl
{

typedef enum
{
  ODOM_MODEL_DIFF,
  ODOM_MODEL_OMNI,
  ODOM_MODEL_DIFF_CORRECTED,
  ODOM_MODEL_OMNI_CORRECTED
} odom_model_t;

// Odometric sensor data
struct AMCLOdomData : public AMCLSensorData
{
  // Odometric pose
  pf_vector_t pose;

  // Change in odometric pose
  pf_vector_t delta;
};

using AMCLOdomDataPtr = std::shared_ptr<AMCLOdomData>;

// Odometric sensor model
class AMCLOdom : public AMCLSensor
{
public:
  // Default constructor
  AMCLOdom();

  void SetModelDiff(double alpha1,
                            double alpha2,
                            double alpha3,
                            double alpha4);

  void SetModelOmni(double alpha1,
                            double alpha2,
                            double alpha3,
                            double alpha4,
                            double alpha5);

  void SetModel( odom_model_t type,
                         double alpha1,
                         double alpha2,
                         double alpha3,
                         double alpha4,
                         double alpha5 = 0 );

  // Update the filter based on the action model.  Returns true if the filter
  // has been updated.
  virtual bool UpdateAction(pf_t *pf, AMCLSensorDataPtr data);

private:
  // Current data timestamp
  double time;

  // Model type
  odom_model_t model_type;

  // Drift parameters
  double alpha1, alpha2, alpha3, alpha4, alpha5;
};

using AMCLOdomPtr = std::shared_ptr<AMCLOdom>;

} // namespace amcl
