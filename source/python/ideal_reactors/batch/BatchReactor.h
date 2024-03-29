#pragma once

// OpenSMOKE++ Definitions
#include <OpenSMOKEpp>

// CHEMKIN maps
#include <maps/Maps_CHEMKIN>

// Utilities
#include <idealreactors/utilities/Utilities>
#include <utilities/ropa/OnTheFlyROPA.h>
#include <utilities/cema/OnTheFlyCEMA.h>
#include <utilities/ontheflypostprocessing/OnTheFlyPostProcessing.h>
#include <utilities/kineticsmodifier/KineticsModifier.h>

// PolimiSoot Analyzer
#include <utilities/soot/polimi/OpenSMOKE_PolimiSoot_Analyzer.h>

// Batch reactor
#include <idealreactors/batch/BatchReactor>

// pyBIND11
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

/**
 * @class BatchReactor
 * @brief HighLevel wrap to the Batch reactor solver of OpenSMOKEpp. This is not a wrap
 * around the Base class but provide the access to most of the functionalities of the
 * BatchReactor class of OpenSMOKEpp.
 */
class BatchReactor {
 public:
  /**
   * @brief Default Constructor. Set some of the default values.
   */
  BatchReactor();

  /**
   * @brief Default destructor.
   */
  ~BatchReactor();

  /**
   * @brief Convenient function that calls the numerical solvers and expose the solution
   * of the reactor.
   */
  void Solve();

  // Setter functions
  /**
   * @brief Function that sets the thermodynamic map onto which performs the simulation.
   * @param thermo Pointer to the thermodynamic map.
   */
  const void SetThermodynamic(OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermo) {
    thermodynamicsMapXML_ = thermo;
  };

  /**
   * @brief Function that sets the kinetics map onto which performs the simulation.
   * @param kinetics Pointer to the kinetic map.
   */
  const void SetKinetics(OpenSMOKE::KineticsMap_CHEMKIN* kinetics) {
    kineticsMapXML_ = kinetics;
  };

  /**
   * @brief Function that activate the storage of the entire time series solution.
   * @param save Boolean value needed for the activation.
   */
  const void SaveResults(const bool& save) { save_ = save; };

  /**
   * @brief Set the initial Temperature of the simulation.
   * @param value Temperature value.
   * @param units Unit of measurements for the Temperature, available are: K | C.
   */
  const void SetTemperature(const double& value, const std::string& units);

  /**
   * @brief Set the initial Pressure of the simulation.
   * @param value Pressure value.
   * @param units Unit of measurements for the Pressure, available are: Pa | bar | atm.
   */
  const void SetPressure(const double& value, const std::string& units);

  /**
   * @brief Set the density of the reacting mixture at the beginning of the simulation.
   * @param value Density value.
   * @param units Unit of measurements for the density, available are: kg/m3 | g/cm3.
   */
  const void SetDensity(const double& value, const std::string& units);

  // Composition
  const void SetInitialComposition(const std::string& initial_composition_type,
                                   const std::vector<std::string>& names,
                                   const std::vector<double>& values);

  // Composition
  const void SetInitialComposition(const std::string& initial_composition_type,
                                   const std::vector<double>& equivalence_ratios,
                                   const std::string& fuel_composition_type,
                                   const std::vector<std::string>& fuel_names,
                                   std::vector<double>& values_fuel,
                                   const std::string& oxidizer_composition_type,
                                   std::vector<std::string>& oxidizer_names,
                                   std::vector<double>& values_oxidizer);

  // End time
  const void SetEndTime(const double& value, const std::string& units);

  // Read start time
  const void SetStartTime(const double& value, const std::string& units);

  // Read volume
  const void SetVolume(const double& value, const std::string& units);

  // Exchange area
  const void SetExchangeArea(const double& value, const std::string& units);

  // Read global thermal exchange coefficient
  const void Set_global_thermal_exchange_coefficient(const double& value,
                                                     const std::string& units);

  // Environment temperature
  const void SetEnvironmentTemperature(const double& value, const std::string& units);

  // Batch reactors Type
  const void SetType(const std::string& value);

  // Set Batch Output options
  const void SetBatchOptions(const bool& verbose = false, const bool& save_results = false,
                             const std::string& output_path = "/dev/null");

  // Set Ode Options
  const void SetOdeOptions(const double& abs_tol, const double& rel_tol);

  // Getter functions
  const double& Tf() const { return Tf_; };

  const double& Pf() const { return Pf_; };

  const Eigen::VectorXd& omegaf() const { return omegaf_; };

  const Eigen::VectorXd& xf() const { return xf_; };

  const std::vector<double>& tOut() const { return tOut_; };

  const std::vector<std::vector<double>>& xOut() const { return xOut_; };

  // Python Wrapper
  const static void BatchReactor_wrapper(py::module_&);

 private:
  const void CleanMemory();

  const void SetAdditionalOptions();

  std::unique_ptr<OpenSMOKE::BatchReactor> batch_;

  OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML_;

  OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML_;

  // Options
  OpenSMOKE::BatchReactor_Options* batch_options_;

  // ODE Parameters
  OpenSMOKE::ODE_Parameters* ode_parameters_;

  // Sensitivity Option
  OpenSMOKE::SensitivityAnalysis_Options* sensitivity_options_;

  // On the fly ROPA
  OpenSMOKE::OnTheFlyROPA* onTheFlyROPA_;

  // On the fly CEMA
  OpenSMOKE::OnTheFlyCEMA* onTheFlyCEMA_;

  // On the fly PostProcessing
  OpenSMOKE::OnTheFlyPostProcessing* on_the_fly_post_processing_;

  // Ignition Delay Times
  OpenSMOKE::IgnitionDelayTimes_Analyzer* idt_;

  OpenSMOKE::BatchReactor_VolumeProfile* batchreactor_volumeprofile_;

  // Polimi soot
  OpenSMOKE::PolimiSoot_Analyzer* polimi_soot_;

  double T, P_Pa;
  OpenSMOKE::OpenSMOKEVectorDouble omega;
  double Tf_, Pf_;
  Eigen::VectorXd omegaf_;
  Eigen::VectorXd xf_;

  std::vector<double> tOut_;
  std::vector<std::vector<double>> xOut_;

  double tEnd_;
  double tStart_;
  double volume_;

  unsigned int state_variables_;
  bool temperature_assigned_;
  bool pressure_assigned_;
  bool density_assigned_;

  double rho_;
  double exchange_area_;
  double global_thermal_exchange_coefficient_;
  double T_environment_;

  bool verbose_;
  bool save_;

  OpenSMOKE::BatchReactor_Type type_;

  const void CeckStatusOfGasMixture() const {
    // TODO refactor this controller.
    if (state_variables_ <= 2) {
      OpenSMOKE::FatalErrorMessage(
          "The status of a gas mixture requires any 2 (and only 2) among: Temperature, "
          "Pressure and Density");
    }
  }
};
