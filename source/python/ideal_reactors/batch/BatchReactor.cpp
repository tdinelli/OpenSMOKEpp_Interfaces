#include "BatchReactor.h"

#include <memory>

#include "math/OpenSMOKEFunctions.h"

BatchReactor::BatchReactor() {
  tStart_ = 0.;  // default 0
  volume_ = 1.;  // default value [1 m3]

  state_variables_ = 0;
  temperature_assigned_ = false;
  pressure_assigned_ = false;
  density_assigned_ = false;

  exchange_area_ = 0.;
  global_thermal_exchange_coefficient_ = 0.;

  verbose_ = true;
};

BatchReactor::~BatchReactor(){};

const void BatchReactor::SetTemperature(const double &value, const std::string &units) {
  if (units == "K") {
    T = value;
  } else if (units == "C") {
    T = value + 273.15;
  } else {
    OpenSMOKE::FatalErrorMessage("Unknown temperature units");
  }
  T_environment_ = T;
  state_variables_++;
  temperature_assigned_ = true;
}

const void BatchReactor::SetPressure(const double &value, const std::string &units) {
  if (units == "Pa") {
    P_Pa = value;
  } else if (units == "bar") {
    P_Pa = value * 1.e5;
  } else if (units == "atm") {
    P_Pa = value * 101325.;
  } else {
    OpenSMOKE::FatalErrorMessage("Unknown pressure units");
  }
  state_variables_++;
  pressure_assigned_ = true;
}

const void BatchReactor::SetDensity(const double &value, const std::string &units) {
  if (units == "kg/m3") {
    rho_ = value;
  } else if (units == "g/cm3") {
    rho_ = value * 1.e3;
  } else {
    OpenSMOKE::FatalErrorMessage("Unknown density units");
  }
  state_variables_++;
  density_assigned_ = true;
}

const void BatchReactor::SetInitialComposition(
    const std::string &initial_composition_type, const std::vector<std::string> &names,
    const std::vector<double> &values) {
  if (initial_composition_type == "MoleFractions") {
    const double sum = std::accumulate(values.begin(), values.end(), 0.);
    if (sum < (1. - 1e-6) || sum > (1. + 1e-6)) {
      OpenSMOKE::FatalErrorMessage("The mole fractions must sum to 1.");
    }

    OpenSMOKE::OpenSMOKEVectorDouble x(thermodynamicsMapXML_->NumberOfSpecies());
    for (unsigned int i = 0; i < names.size(); i++) {
      x[thermodynamicsMapXML_->IndexOfSpecies(names[i])] = values[i] / sum;
    }

    ChangeDimensions(thermodynamicsMapXML_->NumberOfSpecies(), &omega, true);
    double MW;
    thermodynamicsMapXML_->MassFractions_From_MoleFractions(omega.GetHandle(), MW,
                                                            x.GetHandle());
  } else if (initial_composition_type == "MassFractions") {
    const double sum = std::accumulate(values.begin(), values.end(), 0.);
    if (sum < (1. - 1e-6) || sum > (1. + 1e-6)) {
      OpenSMOKE::FatalErrorMessage("The mass fractions must sum to 1.");
    }

    ChangeDimensions(thermodynamicsMapXML_->NumberOfSpecies(), &omega, true);
    for (unsigned int i = 0; i < names.size(); i++) {
      omega[thermodynamicsMapXML_->IndexOfSpecies(names[i])] = values[i] / sum;
    }
  } else if (initial_composition_type == "Moles") {
    const double sum = std::accumulate(values.begin(), values.end(), 0.);
    OpenSMOKE::OpenSMOKEVectorDouble x(thermodynamicsMapXML_->NumberOfSpecies());
    for (unsigned int i = 0; i < names.size(); i++) {
      x[thermodynamicsMapXML_->IndexOfSpecies(names[i])] = values[i] / sum;
    }
    ChangeDimensions(thermodynamicsMapXML_->NumberOfSpecies(), &omega, true);
    double MW;
    thermodynamicsMapXML_->MassFractions_From_MoleFractions(omega.GetHandle(), MW,
                                                            x.GetHandle());
  } else if (initial_composition_type == "Masses") {
    const double sum = std::accumulate(values.begin(), values.end(), 0.);
    ChangeDimensions(thermodynamicsMapXML_->NumberOfSpecies(), &omega, true);
    for (unsigned int i = 0; i < names.size(); i++) {
      omega[thermodynamicsMapXML_->IndexOfSpecies(names[i])] = values[i] / sum;
    }
  }
}

const void BatchReactor::SetInitialComposition(
    const std::string &initial_composition_type,
    const std::vector<double> &equivalence_ratios,
    const std::string &fuel_composition_type, const std::vector<std::string> &fuel_names,
    std::vector<double> &values_fuel, const std::string &oxidizer_composition_type,
    std::vector<std::string> &oxidizer_names, std::vector<double> &values_oxidizer) {
  if (initial_composition_type == "EquivalenceRatio") {
    if (equivalence_ratios.size() != 1) {
      OpenSMOKE::FatalErrorMessage("Only a single EquivalenceRatio can be specified");
    }

    double equivalence_ratio = equivalence_ratios[0];

    if (fuel_composition_type == "FuelMoleFractions") {
      const double sum = std::accumulate(values_fuel.begin(), values_fuel.end(), 0.);
      if (sum < (1. - 1e-6) || sum > (1. + 1e-6)) {
        OpenSMOKE::FatalErrorMessage("The fuel mass fractions must sum to 1.");
      }

      for (unsigned int i = 0; i < fuel_names.size(); i++) {
        values_fuel[i] /= sum;
      }
    } else if (fuel_composition_type == "FuelMassFractions") {
      const double sum = std::accumulate(values_fuel.begin(), values_fuel.end(), 0.);
      if (sum < (1. - 1e-6) || sum > (1. + 1e-6)) {
        OpenSMOKE::FatalErrorMessage("The fuel mass fractions must sum to 1.");
      }

      OpenSMOKE::OpenSMOKEVectorDouble omega_fuel(
          thermodynamicsMapXML_->NumberOfSpecies());
      for (unsigned int i = 0; i < fuel_names.size(); i++) {
        omega_fuel[thermodynamicsMapXML_->IndexOfSpecies(fuel_names[i])] =
            values_fuel[i] / sum;
      }

      double MW_fuel;
      OpenSMOKE::OpenSMOKEVectorDouble x_fuel(thermodynamicsMapXML_->NumberOfSpecies());
      thermodynamicsMapXML_->MoleFractions_From_MassFractions(x_fuel.GetHandle(), MW_fuel,
                                                              omega_fuel.GetHandle());
      for (unsigned int i = 0; i < fuel_names.size(); i++) {
        values_fuel[i] = x_fuel[thermodynamicsMapXML_->IndexOfSpecies(fuel_names[i])];
      }
    } else if (fuel_composition_type == "FuelMoles") {
      // TODO
      // dictionary.ReadOption("@FuelMoles", names_fuel, values_fuel);
    } else if (fuel_composition_type == "FuelMasses") {
      OpenSMOKE::OpenSMOKEVectorDouble omega_fuel(
          thermodynamicsMapXML_->NumberOfSpecies());
      for (unsigned int i = 0; i < fuel_names.size(); i++) {
        omega_fuel[thermodynamicsMapXML_->IndexOfSpecies(fuel_names[i])] = values_fuel[i];
      }

      double MW_fuel;
      OpenSMOKE::OpenSMOKEVectorDouble x_fuel(thermodynamicsMapXML_->NumberOfSpecies());
      thermodynamicsMapXML_->MoleFractions_From_MassFractions(x_fuel.GetHandle(), MW_fuel,
                                                              omega_fuel.GetHandle());
      for (unsigned int i = 0; i < fuel_names.size(); i++) {
        values_fuel[i] = x_fuel[thermodynamicsMapXML_->IndexOfSpecies(fuel_names[i])];
      }
    } else {
      OpenSMOKE::FatalErrorMessage(
          "The EquivalenceRatio option requires the user specifies the fuel "
          "composition");
    }

    if (oxidizer_composition_type == "OxidizerMoleFractions") {
      const double sum =
          std::accumulate(values_oxidizer.begin(), values_oxidizer.end(), 0.);
      if (sum < (1. - 1e-6) || sum > (1. + 1e-6)) {
        OpenSMOKE::FatalErrorMessage("The oxidizer mass fractions must sum to 1.");
      }

      for (unsigned int i = 0; i < oxidizer_names.size(); i++) {
        values_oxidizer[i] /= sum;
      }
    } else if (oxidizer_composition_type == "OxidizerMassFractions") {
      const double sum =
          std::accumulate(values_oxidizer.begin(), values_oxidizer.end(), 0.);
      if (sum < (1. - 1e-6) || sum > (1. + 1e-6)) {
        OpenSMOKE::FatalErrorMessage("The oxidizer mass fractions must sum to 1.");
      }

      OpenSMOKE::OpenSMOKEVectorDouble omega_oxidizer(
          thermodynamicsMapXML_->NumberOfSpecies());
      for (unsigned int i = 0; i < oxidizer_names.size(); i++) {
        omega_oxidizer[thermodynamicsMapXML_->IndexOfSpecies(oxidizer_names[i])] =
            values_oxidizer[i] / sum;
      }

      double MW_oxidizer;
      OpenSMOKE::OpenSMOKEVectorDouble x_oxidizer(
          thermodynamicsMapXML_->NumberOfSpecies());
      thermodynamicsMapXML_->MoleFractions_From_MassFractions(
          x_oxidizer.GetHandle(), MW_oxidizer, omega_oxidizer.GetHandle());
      for (unsigned int i = 0; i < oxidizer_names.size(); i++) {
        values_oxidizer[i] =
            x_oxidizer[thermodynamicsMapXML_->IndexOfSpecies(oxidizer_names[i])];
      }
    } else if (oxidizer_composition_type == "OxidizerMoles") {
      // TODO
      // dictionary.ReadOption("@OxidizerMoles", names_oxidizer, values_oxidizer);
    } else if (oxidizer_composition_type == "OxidizerMasses") {
      OpenSMOKE::OpenSMOKEVectorDouble omega_oxidizer(
          thermodynamicsMapXML_->NumberOfSpecies());
      for (unsigned int i = 0; i < oxidizer_names.size(); i++) {
        omega_oxidizer[thermodynamicsMapXML_->IndexOfSpecies(oxidizer_names[i])] =
            values_oxidizer[i];
      }

      double MW_oxidizer;
      OpenSMOKE::OpenSMOKEVectorDouble x_oxidizer(
          thermodynamicsMapXML_->NumberOfSpecies());
      thermodynamicsMapXML_->MoleFractions_From_MassFractions(
          x_oxidizer.GetHandle(), MW_oxidizer, omega_oxidizer.GetHandle());
      for (unsigned int i = 0; i < oxidizer_names.size(); i++) {
        values_oxidizer[i] =
            x_oxidizer[thermodynamicsMapXML_->IndexOfSpecies(oxidizer_names[i])];
      }
    } else {
      oxidizer_names.resize(2);
      oxidizer_names[0] = "O2";
      oxidizer_names[1] = "N2";
      values_oxidizer.resize(2);
      values_oxidizer[0] = 0.21;
      values_oxidizer[1] = 0.79;
    }

    std::vector<double> values =
        thermodynamicsMapXML_->GetMoleFractionsFromEquivalenceRatio(
            equivalence_ratio, fuel_names, values_fuel, oxidizer_names, values_oxidizer);

    const double sum = std::accumulate(values.begin(), values.end(), 0.);
    if (sum < (1. - 1e-6) || sum > (1. + 1e-6)) {
      OpenSMOKE::FatalErrorMessage("The mole fractions must sum to 1.");
    }

    OpenSMOKE::OpenSMOKEVectorDouble x(thermodynamicsMapXML_->NumberOfSpecies());
    for (unsigned int i = 0; i < thermodynamicsMapXML_->NumberOfSpecies(); i++) {
      x[i + 1] = values[i] / sum;
    }
    ChangeDimensions(thermodynamicsMapXML_->NumberOfSpecies(), &omega, true);
    double MW;
    thermodynamicsMapXML_->MassFractions_From_MoleFractions(omega.GetHandle(), MW,
                                                            x.GetHandle());
  }

  if (density_assigned_ == true) {
    const double MW =
        thermodynamicsMapXML_->MolecularWeight_From_MassFractions(omega.GetHandle());
    if (temperature_assigned_ == true) {
      P_Pa = rho_ * PhysicalConstants::R_J_kmol * T / MW;
    }
    if (pressure_assigned_ == true) {
      T = P_Pa * MW / PhysicalConstants::R_J_kmol / rho_;
    }
  }
}

const void BatchReactor::SetEndTime(const double &value, const std::string &units) {
  if (units == "s") {
    tEnd_ = value;
  } else if (units == "ms") {
    tEnd_ = value / 1000.;
  } else if (units == "min") {
    tEnd_ = value * 60.;
  } else if (units == "h") {
    tEnd_ = value * 3600.;
  } else {
    OpenSMOKE::FatalErrorMessage("Unknown time units");
  }
}

const void BatchReactor::SetStartTime(const double &value, const std::string &units) {
  if (units == "s") {
    tStart_ = value;
  } else if (units == "ms") {
    tStart_ = value / 1000.;
  } else if (units == "min") {
    tStart_ = value * 60.;
  } else if (units == "h") {
    tStart_ = value * 3600.;
  } else {
    OpenSMOKE::FatalErrorMessage("Unknown time units");
  }
}

const void BatchReactor::SetVolume(const double &value, const std::string &units) {
  if (units == "m3") {
    volume_ = value;
  } else if (units == "dm3") {
    volume_ = value / 1.e3;
  } else if (units == "cm3") {
    volume_ = value / 1.e6;
  } else if (units == "mm3") {
    volume_ = value / 1.e9;
  } else if (units == "l") {
    volume_ = value / 1.e3;
  } else {
    OpenSMOKE::FatalErrorMessage("Unknown volume_ units");
  }
}

const void BatchReactor::SetExchangeArea(const double &value, const std::string &units) {
  if (units == "m2") {
    exchange_area_ = value;
  } else if (units == "dm2") {
    exchange_area_ = value / 1.e2;
  } else if (units == "cm2") {
    exchange_area_ = value / 1.e4;
  } else if (units == "mm2") {
    exchange_area_ = value / 1.e6;
  } else {
    OpenSMOKE::FatalErrorMessage("Unknown area units");
  }
}

const void BatchReactor::Set_global_thermal_exchange_coefficient(
    const double &value, const std::string &units) {
  if (units == "W/m2/K") {
    global_thermal_exchange_coefficient_ = value;
  } else if (units == "W/m2/C") {
    global_thermal_exchange_coefficient_ = value;
  } else if (units == "kcal/m2/K") {
    global_thermal_exchange_coefficient_ = value * 4186.8;
  } else if (units == "kcal/m2/C") {
    global_thermal_exchange_coefficient_ = value * 4186.8;
  } else {
    OpenSMOKE::FatalErrorMessage("Unknown global thermal exchange coefficient units");
  }
}

const void BatchReactor::SetEnvironmentTemperature(const double &value,
                                                   const std::string &units) {
  if (units == "K") {
    T_environment_ = value;
  } else if (units == "C") {
    T_environment_ = value + 273.15;
  } else {
    OpenSMOKE::FatalErrorMessage("Unknown temperature units");
  }
}

const void BatchReactor::SetType(const std::string &value) {
  if (value == "Isothermal-ConstantVolume") {
    type_ = OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTV;
  } else if (value == "Isothermal-ConstantPressure") {
    type_ = OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTP;
  } else if (value == "NonIsothermal-ConstantVolume") {
    type_ = OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_CONSTANTV;
  } else if (value == "NonIsothermal-ConstantPressure") {
    type_ = OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_CONSTANTP;
  } else if (value == "NonIsothermal-UserDefinedVolume") {
    type_ = OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_USERDEFINEDVOLUME;
  } else {
    OpenSMOKE::FatalErrorMessage("Unknown batch reactor type: " + value);
  }
}

const void BatchReactor::SetBatchOptions(const bool &verbose, const bool &save_results,
                                         const std::string &output_path) {
  batch_options_ = new OpenSMOKE::BatchReactor_Options();

  const boost::filesystem::path dummy_path = output_path;

  batch_options_->SetVerboseVideo(verbose);
  batch_options_->SetOutputPath(dummy_path);
  batch_options_->SetNumberOfSteps_File(1);
}

const void BatchReactor::SetOdeOptions(const double &abs_tol, const double &rel_tol) {
  ode_parameters_ = new OpenSMOKE::ODE_Parameters();

  ode_parameters_->SetRelativeTolerance(abs_tol);
  ode_parameters_->SetAbsoluteTolerance(rel_tol);
  // ode_parameters_->SetMinimumStep();
  // ode_parameters_->SetMaximumStep();
  // ode_parameters_->SetInitialStep();
  // ode_parameters_->SetMaximumNumberOfSteps();
  // ode_parameters_->SetMaximumOrder();
  // ode_parameters_->SetFullPivoting();
}

const void BatchReactor::SetAdditionalOptions() {
  sensitivity_options_ = new OpenSMOKE::SensitivityAnalysis_Options();

  onTheFlyROPA_ = new OpenSMOKE::OnTheFlyROPA(*thermodynamicsMapXML_, *kineticsMapXML_);

  onTheFlyCEMA_ = new OpenSMOKE::OnTheFlyCEMA(*thermodynamicsMapXML_, *kineticsMapXML_,
                                              batch_options_->output_path());

  on_the_fly_post_processing_ = new OpenSMOKE::OnTheFlyPostProcessing(
      *thermodynamicsMapXML_, *kineticsMapXML_, batch_options_->output_path());

  idt_ = new OpenSMOKE::IgnitionDelayTimes_Analyzer();

  batchreactor_volumeprofile_ = new OpenSMOKE::BatchReactor_VolumeProfile();

  polimi_soot_ = new OpenSMOKE::PolimiSoot_Analyzer(thermodynamicsMapXML_);
}

const void BatchReactor::CleanMemory() {
  // Cleaning state variables in order to handle easily multiple serial simulations
  // state_variables_ = 0;

  // At the moment clean memory handles only the additional not modified/used stuff
  delete sensitivity_options_;
  sensitivity_options_ = NULL;

  delete onTheFlyROPA_;
  onTheFlyROPA_ = NULL;

  delete onTheFlyCEMA_;
  onTheFlyCEMA_ = NULL;

  delete on_the_fly_post_processing_;
  on_the_fly_post_processing_ = NULL;

  delete idt_;

  delete batchreactor_volumeprofile_;
  batchreactor_volumeprofile_ = NULL;

  delete polimi_soot_;
  polimi_soot_ = NULL;
}

void BatchReactor::Solve() {
  // Perform check of the thermochemistry state of the system
  CeckStatusOfGasMixture();

  // Solve the ODE system: NonIsothermal, Constant Volume
  if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_CONSTANTV) {
    batch_ = std::make_unique<OpenSMOKE::BatchReactor_NonIsothermal_ConstantVolume>(
        *thermodynamicsMapXML_, *kineticsMapXML_, *ode_parameters_, *batch_options_,
        *onTheFlyROPA_, *onTheFlyCEMA_, *on_the_fly_post_processing_, *idt_,
        *polimi_soot_, volume_, T, P_Pa, omega, global_thermal_exchange_coefficient_,
        exchange_area_, T_environment_);

    batch_->SaveResults(save_);
    batch_->Solve(tStart_, tEnd_);
    tOut_ = batch_->tOut();
    xOut_ = batch_->xOut();
  } else if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_USERDEFINEDVOLUME) {
    std::cout << "USER DEFINED VOLUME batch reactor not implemented yet!" << std::endl;
    exit(-1);
  } else if (type_ == OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTV) {
    batch_ = std::make_unique<OpenSMOKE::BatchReactor_Isothermal_ConstantVolume>(
        *thermodynamicsMapXML_, *kineticsMapXML_, *ode_parameters_, *batch_options_,
        *onTheFlyROPA_, *onTheFlyCEMA_, *on_the_fly_post_processing_, *idt_,
        *polimi_soot_, volume_, T, P_Pa, omega);

    batch_->Solve(tStart_, tEnd_);
  } else if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_CONSTANTP) {
    batch_ = std::make_unique<OpenSMOKE::BatchReactor_NonIsothermal_ConstantPressure>(
        *thermodynamicsMapXML_, *kineticsMapXML_, *ode_parameters_, *batch_options_,
        *onTheFlyROPA_, *onTheFlyCEMA_, *on_the_fly_post_processing_, *idt_,
        *polimi_soot_, volume_, T, P_Pa, omega, global_thermal_exchange_coefficient_,
        exchange_area_, T_environment_);
    batch_->Solve(tStart_, tEnd_);
  } else if (type_ == OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTP) {
    batch_ = std::make_unique<OpenSMOKE::BatchReactor_Isothermal_ConstantPressure>(
        *thermodynamicsMapXML_, *kineticsMapXML_, *ode_parameters_, *batch_options_,
        *onTheFlyROPA_, *onTheFlyCEMA_, *on_the_fly_post_processing_, *idt_,
        *polimi_soot_, volume_, T, P_Pa, omega);
    batch_->Solve(tStart_, tEnd_);
  } else {
    OpenSMOKE::FatalErrorMessage("Unknown batch reactor type or type not setted");
  }

  OpenSMOKE::OpenSMOKEVectorDouble tmp_omegaf_(thermodynamicsMapXML_->NumberOfSpecies());
  OpenSMOKE::OpenSMOKEVectorDouble tmp_xf_(thermodynamicsMapXML_->NumberOfSpecies());

  batch_->GetFinalStatus(Tf_, Pf_, tmp_omegaf_);

  double MWf_ =
      thermodynamicsMapXML_->MolecularWeight_From_MassFractions(tmp_omegaf_.GetHandle());

  thermodynamicsMapXML_->MoleFractions_From_MassFractions(tmp_xf_.GetHandle(), MWf_,
                                                          tmp_omegaf_.GetHandle());

  xf_.resize(tmp_xf_.Size());
  tmp_xf_.CopyTo(xf_.data());

  omegaf_.resize(tmp_omegaf_.Size());
  tmp_omegaf_.CopyTo(omegaf_.data());
}

const void BatchReactor::BatchReactor_wrapper(py::module_ &m) {
  constexpr auto call_guard = py::call_guard<py::gil_scoped_release>();

  py::class_<BatchReactor>(m, "BatchReactor")
      .def(py::init<>(), call_guard, "Constructor")
      .def("CleanMemory", &BatchReactor::CleanMemory, call_guard, "")
      .def("SetThermodynamic", &BatchReactor::SetThermodynamic, call_guard,
           "Set thermodynamic map object")
      .def("SetKinetics", &BatchReactor::SetKinetics, call_guard, "")
      .def("SetTemperature", &BatchReactor::SetTemperature, call_guard, "")
      .def("SaveResults", &BatchReactor::SaveResults, call_guard, "")
      .def("SetPressure", &BatchReactor::SetPressure, call_guard, "")
      .def("SetDensity", &BatchReactor::SetDensity, call_guard, "")
      .def("SetStartTime", &BatchReactor::SetStartTime, call_guard, "")
      .def("SetEndTime", &BatchReactor::SetEndTime, call_guard, "")
      .def("SetVolume", &BatchReactor::SetVolume, call_guard, "")
      .def("SetExchangeArea", &BatchReactor::SetExchangeArea, call_guard, "")
      .def("Set_global_thermal_exchange_coefficient",
           &BatchReactor::Set_global_thermal_exchange_coefficient, call_guard, "")
      .def("SetEnvironmentTemperature", &BatchReactor::SetEnvironmentTemperature,
           call_guard, "")
      .def("SetType", &BatchReactor::SetType, call_guard, "")
      .def("SetBatchOptions", &BatchReactor::SetBatchOptions, call_guard, "",
           py::arg("verbose") = false, py::arg("save_results") = false,
           py::arg("output_path") = "/dev/null")
      .def("SetOdeOptions", &BatchReactor::SetOdeOptions, call_guard, "")
      .def("SetAdditionalOptions", &BatchReactor::SetAdditionalOptions, call_guard, "")
      .def("SetInitialComposition",
           py::overload_cast<const std::string &, const std::vector<std::string> &,
                             const std::vector<double> &>(
               &BatchReactor::SetInitialComposition),
           call_guard, "")
      .def("SetInitialComposition",
           py::overload_cast<const std::string &, const std::vector<double> &,
                             const std::string &, const std::vector<std::string> &,
                             std::vector<double> &, const std::string &,
                             std::vector<std::string> &, std::vector<double> &>(
               &BatchReactor::SetInitialComposition),
           call_guard, "")
      .def("Tf", &BatchReactor::Tf, call_guard, "")
      .def("Pf", &BatchReactor::Pf, call_guard, "")
      .def("omegaf", &BatchReactor::omegaf, call_guard, "")
      .def("xf", &BatchReactor::xf, call_guard, "")
      .def("tOut", &BatchReactor::tOut, call_guard, "")
      .def("xOut", &BatchReactor::xOut, call_guard, "")
      .def("Solve", &BatchReactor::Solve, call_guard, "");
}
