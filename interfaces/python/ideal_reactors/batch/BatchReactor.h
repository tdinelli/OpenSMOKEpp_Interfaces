#pragma once

// clang-format off
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
#include "maps/KineticsMap_CHEMKIN.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
// clang-format on

namespace py = pybind11;

class BatchReactor {
 public:
  BatchReactor();

  ~BatchReactor();

  void Solve();

  // Setter functions
  const void SetThermodynamic(OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermo) {
    thermodynamicsMapXML_ = thermo;
  };

  const void SetKinetics(OpenSMOKE::KineticsMap_CHEMKIN* kinetics) {
    kineticsMapXML_ = kinetics;
  };

  // Temperature
  const void SetTemperature(const double& value, const std::string& units);

  // Pressure
  const void SetPressure(const double& value, const std::string& units);

  // Density
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
  const void SetBatchOptions();

  // Set Ode Options
  const void SetOdeOptions();
  // const std::vector<double>& GetTime() const { return time_vector_; };
  // const std::vector<double>& GetTemperature() const { return temperature_vector_; };
  // const std::vector<double>& GetPressure() const { return pressure_vector_; };
  //
  // const std::vector<std::vector<double>>& GetMoleFractions(const
  // std::vector<std::string> &species_names) {
  //   std::vector<std::vector<double>> selected_molefractions(species_names.size(),
  //   std::vector<double>(time_vector_.size())); for (unsigned int i = 0; i <
  //   species_names.size(); i++) {
  //     unsigned int j = thermodynamicsMapXML->IndexOfSpecies(species_names[i]);
  //     for (unsigned int k = 0; k < time_vector_.size(); k++)
  //       selected_molefractions[i][k] = mole_fractions_[k][j];
  //   }
  //   return selected_molefractions;
  // }

  // std::vector<std::vector<double>> GetMassFractions(
  //     std::vector<std::string> species_names) {
  //   std::vector<std::vector<double>> selected_massfractions(
  //       species_names.size(), std::vector<double>(time_vector_.size()));
  //   for (unsigned int i = 0; i < species_names.size(); i++) {
  //     unsigned int j = thermodynamicsMapXML->IndexOfSpecies(species_names[i]);
  //     for (unsigned int k = 0; k < time_vector_.size(); k++)
  //       selected_massfractions[i][k] = mass_fractions_[k][j];
  //   }
  //   return selected_massfractions;
  // }

  // Getter functions

  // Python Wrapper
  const static void BatchReactor_wrapper(py::module_&);

 private:
  OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML_;

  OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML_;

  // Options
  OpenSMOKE::BatchReactor_Options* batch_options_;

  // ODE Parameters
  OpenSMOKE::ODE_Parameters* ode_parameters_;

  // Sensitivity Options
  OpenSMOKE::SensitivityAnalysis_Options* sensitivity_options_;

  // On the fly ROPA
  OpenSMOKE::OnTheFlyROPA*
      onTheFlyROPA_;  // new OpenSMOKE::OnTheFlyROPA(*thermodynamicsMapXML,
                      // *kineticsMapXML);

  // On the fly CEMA
  OpenSMOKE::OnTheFlyCEMA*
      onTheFlyCEMA_;  // = new OpenSMOKE::OnTheFlyCEMA( *thermodynamicsMapXML,
                      // *kineticsMapXML, batch_options->output_path());

  // On the fly PostProcessing
  OpenSMOKE::OnTheFlyPostProcessing*
      on_the_fly_post_processing_;  // = new
                                    // OpenSMOKE::OnTheFlyPostProcessing(*thermodynamicsMapXML,
                                    // *kineticsMapXML, batch_options->output_path());

  // Ignition Delay Times
  OpenSMOKE::IgnitionDelayTimes_Analyzer*
      idt_;  // = new OpenSMOKE::IgnitionDelayTimes_Analyzer();

  OpenSMOKE::BatchReactor_VolumeProfile* batchreactor_volumeprofile_;

  // Polimi soot
  OpenSMOKE::PolimiSoot_Analyzer*
      polimi_soot_;  // = new OpenSMOKE::PolimiSoot_Analyzer(thermodynamicsMapXML);

  double T, P_Pa;
  double Tf_, Pf_;
  OpenSMOKE::OpenSMOKEVectorDouble omega;
  Eigen::VectorXd omegaf_;
  Eigen::VectorXd xf_;

  double tEnd_;
  double tStart_ = 0.;  // default 0
  double volume_ = 1.;  // default value [1 m3]

  unsigned int state_variables_ = 0;
  bool temperature_assigned_ = false;
  bool pressure_assigned_ = false;
  bool density_assigned_ = false;

  double rho_;
  double exchange_area_ = 0.;
  double global_thermal_exchange_coefficient_ = 0.;
  double T_environment_ = T;

  bool verbose_ = true;
  std::vector<double> time_vector_;
  std::vector<double> temperature_vector_;
  std::vector<double> pressure_vector_;
  std::vector<OpenSMOKE::OpenSMOKEVectorDouble> mole_fractions_;
  std::vector<OpenSMOKE::OpenSMOKEVectorDouble> mass_fractions_;

  OpenSMOKE::BatchReactor_Type type_;

  const void CeckStatusOfGasMixture() {
    std::cout << "State variables: " << state_variables_ << std::endl;
    if (state_variables_ != 2) {
      OpenSMOKE::FatalErrorMessage(
          "The status of a gas mixture requires any 2 (and only 2) among: Temperature, "
          "Pressure and Density");
    }
  }
};
