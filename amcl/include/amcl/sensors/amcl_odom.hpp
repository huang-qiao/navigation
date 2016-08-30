#pragma once

#include "amcl_sensor.hpp"
#include "../pf/pf_pdf.hpp"

namespace amcl
{

enum class OdomModel
{
  DIFF,
  OMNI,
  DIFF_CORRECTED,
  OMNI_CORRECTED
};

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

  void SetModel( OdomModel type,
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
  OdomModel model_type;

  // Drift parameters
  double alpha1, alpha2, alpha3, alpha4, alpha5;
};

using AMCLOdomPtr = std::shared_ptr<AMCLOdom>;

} // namespace amcl
