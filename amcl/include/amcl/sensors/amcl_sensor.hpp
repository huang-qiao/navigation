#pragma once

#include "../pf/pf.hpp"
#include <memory>

namespace amcl
{

// Forward declarations
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
  virtual bool UpdateAction(pf_t *pf, AMCLSensorDataPtr data);

  // Initialize the filter based on the sensor model.  Returns true if the
  // filter has been initialized.
  virtual bool InitSensor(pf_t *pf, AMCLSensorDataPtr data);

  // Update the filter based on the sensor model.  Returns true if the
  // filter has been updated.
  virtual bool UpdateSensor(pf_t *pf, AMCLSensorDataPtr data);

  // Flag is true if this is the action sensor
  bool is_action;

  // Action pose (action sensors only)
  Pose pose;

  // AMCL Base
  //protected: AdaptiveMCL & AMCL;

#ifdef INCLUDE_RTKGUI
  // Setup the GUI
  virtual void SetupGUI(rtk_canvas_t *canvas, rtk_fig_t *robot_fig);

  // Finalize the GUI
  virtual void ShutdownGUI(rtk_canvas_t *canvas, rtk_fig_t *robot_fig);

  // Draw sensor data
  virtual void UpdateGUI(rtk_canvas_t *canvas, rtk_fig_t *robot_fig, AMCLSensorDataPtr data);
#endif
};

using AMCLSensorPtr = std::shared_ptr<AMCLSensor>;

// Base class for all AMCL sensor measurements
struct AMCLSensorData
{
  // Pointer to sensor that generated the data
  AMCLSensorPtr sensor;
  virtual ~AMCLSensorData() {}

  // Data timestamp
  double time;
};

} // namespace amcl;
