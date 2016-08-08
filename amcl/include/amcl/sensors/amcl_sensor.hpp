#pragma once

#include "../pf/pf.hpp"

namespace amcl {

// Forward declarations
class AMCLSensor;
using AMCLSensorPtr = std::shared_ptr<AMCLSensor>;

struct AMCLSensorData;
using AMCLSensorDataPtr = std::shared_ptr<AMCLSensorData>;

// Base class for all AMCL sensors
class AMCLSensor
{
public:
  // Default constructor
  AMCLSensor();
         
  // Default destructor
  virtual ~AMCLSensor();

  // Update the filter based on the action model.  Returns true if the filter
  // has been updated.
  virtual bool updateAction(pf::ParticleFilterPtr pf, AMCLSensorDataPtr data);

  // Initialize the filter based on the sensor model.  Returns true if the
  // filter has been initialized.
  virtual bool initSensor(pf::ParticleFilterPtr pf, AMCLSensorDataPtr data);

  // Update the filter based on the sensor model.  Returns true if the
  // filter has been updated.
  virtual bool updateSensor(pf::ParticleFilterPtr pf, AMCLSensorDataPtr data);

  // Flag is true if this is the action sensor
  bool is_action_;

  // Action pose (action sensors only)
  pf::Pose pose_;
};

// Base class for all AMCL sensor measurements
struct AMCLSensorData
{
public:
  // Pointer to sensor that generated the data
  AMCLSensorPtr sensor;
  virtual ~AMCLSensorData() {}

  // Data timestamp
  double time;
};

} // namespace amcl
