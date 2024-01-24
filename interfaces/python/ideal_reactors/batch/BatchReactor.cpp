#include "BatchReactor.h"

void BatchReactor::SetTemperature(double value, std::string units) {
  if (units == "K") {
    T = value;
  } else if (units == "C") {
    T = value + 273.15;
  } else {
    OpenSMOKE::FatalErrorMessage("Unknown temperature units");
  }
  state_variables_++;
  temperature_assigned_ = true;
}

void BatchReactor::SetPressure(double value, std::string units) {
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

void BatchReactor::SetDensity(double value, std::string units) {
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

void BatchReactor::SetInitialComposition(std::string initial_composition_type,
                                         std::vector<std::string> names,
                                         std::vector<double> values) {
  if (initial_composition_type == "MoleFractions") {
    const double sum = std::accumulate(values.begin(), values.end(), 0.);
    if (sum < (1. - 1e-6) || sum > (1. + 1e-6)) {
      OpenSMOKE::FatalErrorMessage("The mole fractions must sum to 1.");
    }

    OpenSMOKE::OpenSMOKEVectorDouble x(thermodynamicsMapXML->NumberOfSpecies());
    for (unsigned int i = 0; i < names.size(); i++) {
      x[thermodynamicsMapXML->IndexOfSpecies(names[i])] = values[i] / sum;
    }

    ChangeDimensions(thermodynamicsMapXML->NumberOfSpecies(), &omega, true);
    double MW;
    thermodynamicsMapXML->MassFractions_From_MoleFractions(omega.GetHandle(), MW,
                                                           x.GetHandle());
  } else if (initial_composition_type == "MassFractions") {
    const double sum = std::accumulate(values.begin(), values.end(), 0.);
    if (sum < (1. - 1e-6) || sum > (1. + 1e-6)) {
      OpenSMOKE::FatalErrorMessage("The mass fractions must sum to 1.");
    }

    ChangeDimensions(thermodynamicsMapXML->NumberOfSpecies(), &omega, true);
    for (unsigned int i = 0; i < names.size(); i++) {
      omega[thermodynamicsMapXML->IndexOfSpecies(names[i])] = values[i] / sum;
    }
  } else if (initial_composition_type == "Moles") {
    const double sum = std::accumulate(values.begin(), values.end(), 0.);
    OpenSMOKE::OpenSMOKEVectorDouble x(thermodynamicsMapXML->NumberOfSpecies());
    for (unsigned int i = 0; i < names.size(); i++) {
      x[thermodynamicsMapXML->IndexOfSpecies(names[i])] = values[i] / sum;
    }
    ChangeDimensions(thermodynamicsMapXML->NumberOfSpecies(), &omega, true);
    double MW;
    thermodynamicsMapXML->MassFractions_From_MoleFractions(omega.GetHandle(), MW,
                                                           x.GetHandle());
  } else if (initial_composition_type == "Masses") {
    const double sum = std::accumulate(values.begin(), values.end(), 0.);
    ChangeDimensions(thermodynamicsMapXML->NumberOfSpecies(), &omega, true);
    for (unsigned int i = 0; i < names.size(); i++) {
      omega[thermodynamicsMapXML->IndexOfSpecies(names[i])] = values[i] / sum;
    }
  }
}

void BatchReactor::SetInitialComposition(
    std::string initial_composition_type, std::vector<double> equivalence_ratios,
    std::string fuel_composition_type, std::vector<std::string> names_fuel,
    std::vector<double> values_fuel, std::string oxidizer_composition_type,
    std::vector<std::string> names_oxidizer, std::vector<double> values_oxidizer) {
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

      for (unsigned int i = 0; i < names_fuel.size(); i++) {
        values_fuel[i] /= sum;
      }
    } else if (fuel_composition_type == "FuelMassFractions") {
      const double sum = std::accumulate(values_fuel.begin(), values_fuel.end(), 0.);
      if (sum < (1. - 1e-6) || sum > (1. + 1e-6)) {
        OpenSMOKE::FatalErrorMessage("The fuel mass fractions must sum to 1.");
      }

      OpenSMOKE::OpenSMOKEVectorDouble omega_fuel(
          thermodynamicsMapXML->NumberOfSpecies());
      for (unsigned int i = 0; i < names_fuel.size(); i++) {
        omega_fuel[thermodynamicsMapXML->IndexOfSpecies(names_fuel[i])] =
            values_fuel[i] / sum;
      }

      double MW_fuel;
      OpenSMOKE::OpenSMOKEVectorDouble x_fuel(thermodynamicsMapXML->NumberOfSpecies());
      thermodynamicsMapXML->MoleFractions_From_MassFractions(x_fuel.GetHandle(), MW_fuel,
                                                             omega_fuel.GetHandle());
      for (unsigned int i = 0; i < names_fuel.size(); i++) {
        values_fuel[i] = x_fuel[thermodynamicsMapXML->IndexOfSpecies(names_fuel[i])];
      }
    } else if (fuel_composition_type == "FuelMoles") {
      // TODO
      // dictionary.ReadOption("@FuelMoles", names_fuel, values_fuel);
    } else if (fuel_composition_type == "FuelMasses") {
      OpenSMOKE::OpenSMOKEVectorDouble omega_fuel(
          thermodynamicsMapXML->NumberOfSpecies());
      for (unsigned int i = 0; i < names_fuel.size(); i++) {
        omega_fuel[thermodynamicsMapXML->IndexOfSpecies(names_fuel[i])] = values_fuel[i];
      }

      double MW_fuel;
      OpenSMOKE::OpenSMOKEVectorDouble x_fuel(thermodynamicsMapXML->NumberOfSpecies());
      thermodynamicsMapXML->MoleFractions_From_MassFractions(x_fuel.GetHandle(), MW_fuel,
                                                             omega_fuel.GetHandle());
      for (unsigned int i = 0; i < names_fuel.size(); i++) {
        values_fuel[i] = x_fuel[thermodynamicsMapXML->IndexOfSpecies(names_fuel[i])];
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

      for (unsigned int i = 0; i < names_oxidizer.size(); i++) {
        values_oxidizer[i] /= sum;
      }
    } else if (oxidizer_composition_type == "OxidizerMassFractions") {
      const double sum =
          std::accumulate(values_oxidizer.begin(), values_oxidizer.end(), 0.);
      if (sum < (1. - 1e-6) || sum > (1. + 1e-6)) {
        OpenSMOKE::FatalErrorMessage("The oxidizer mass fractions must sum to 1.");
      }

      OpenSMOKE::OpenSMOKEVectorDouble omega_oxidizer(
          thermodynamicsMapXML->NumberOfSpecies());
      for (unsigned int i = 0; i < names_oxidizer.size(); i++) {
        omega_oxidizer[thermodynamicsMapXML->IndexOfSpecies(names_oxidizer[i])] =
            values_oxidizer[i] / sum;
      }

      double MW_oxidizer;
      OpenSMOKE::OpenSMOKEVectorDouble x_oxidizer(
          thermodynamicsMapXML->NumberOfSpecies());
      thermodynamicsMapXML->MoleFractions_From_MassFractions(
          x_oxidizer.GetHandle(), MW_oxidizer, omega_oxidizer.GetHandle());
      for (unsigned int i = 0; i < names_oxidizer.size(); i++) {
        values_oxidizer[i] =
            x_oxidizer[thermodynamicsMapXML->IndexOfSpecies(names_oxidizer[i])];
      }
    } else if (oxidizer_composition_type == "OxidizerMoles") {
      // TODO
      // dictionary.ReadOption("@OxidizerMoles", names_oxidizer, values_oxidizer);
    } else if (oxidizer_composition_type == "OxidizerMasses") {
      OpenSMOKE::OpenSMOKEVectorDouble omega_oxidizer(
          thermodynamicsMapXML->NumberOfSpecies());
      for (unsigned int i = 0; i < names_oxidizer.size(); i++) {
        omega_oxidizer[thermodynamicsMapXML->IndexOfSpecies(names_oxidizer[i])] =
            values_oxidizer[i];
      }

      double MW_oxidizer;
      OpenSMOKE::OpenSMOKEVectorDouble x_oxidizer(
          thermodynamicsMapXML->NumberOfSpecies());
      thermodynamicsMapXML->MoleFractions_From_MassFractions(
          x_oxidizer.GetHandle(), MW_oxidizer, omega_oxidizer.GetHandle());
      for (unsigned int i = 0; i < names_oxidizer.size(); i++) {
        values_oxidizer[i] =
            x_oxidizer[thermodynamicsMapXML->IndexOfSpecies(names_oxidizer[i])];
      }
    } else {
      names_oxidizer.resize(2);
      names_oxidizer[0] = "O2";
      names_oxidizer[1] = "N2";
      values_oxidizer.resize(2);
      values_oxidizer[0] = 0.21;
      values_oxidizer[1] = 0.79;
    }

    std::vector<double> values =
        thermodynamicsMapXML->GetMoleFractionsFromEquivalenceRatio(
            equivalence_ratio, names_fuel, values_fuel, names_oxidizer, values_oxidizer);

    const double sum = std::accumulate(values.begin(), values.end(), 0.);
    if (sum < (1. - 1e-6) || sum > (1. + 1e-6)) {
      OpenSMOKE::FatalErrorMessage("The mole fractions must sum to 1.");
    }

    OpenSMOKE::OpenSMOKEVectorDouble x(thermodynamicsMapXML->NumberOfSpecies());
    for (unsigned int i = 0; i < thermodynamicsMapXML->NumberOfSpecies(); i++) {
      x[i + 1] = values[i] / sum;
    }
    ChangeDimensions(thermodynamicsMapXML->NumberOfSpecies(), &omega, true);
    double MW;
    thermodynamicsMapXML->MassFractions_From_MoleFractions(omega.GetHandle(), MW,
                                                           x.GetHandle());
  }

  if (density_assigned_ == true) {
    const double MW =
        thermodynamicsMapXML->MolecularWeight_From_MassFractions(omega.GetHandle());
    if (temperature_assigned_ == true) {
      P_Pa = rho_ * PhysicalConstants::R_J_kmol * T / MW;
    }
    if (pressure_assigned_ == true) {
      T = P_Pa * MW / PhysicalConstants::R_J_kmol / rho_;
    }
  }
}

void BatchReactor::SetEndTime(double value, std::string units) {
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

void BatchReactor::SetStartTime(double value, std::string units) {
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

void BatchReactor::SetVolume(double value, std::string units) {
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

void BatchReactor::SetExchangeArea(double value, std::string units) {
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

void BatchReactor::Set_global_thermal_exchange_coefficient(double value,
                                                           std::string units) {
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

void BatchReactor::SetEnvironmentTemperature(double value, std::string units) {
  if (units == "K") {
    T_environment_ = value;
  } else if (units == "C") {
    T_environment_ = value + 273.15;
  } else {
    OpenSMOKE::FatalErrorMessage("Unknown temperature units");
  }
}

void BatchReactor::SetType(std::string value) {
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
    OpenSMOKE::FatalErrorMessage("Unknown batch reactor type_: " + value);
  }
}

void BatchReactor::SetBatchOptions(std::string output_path) {
  batch_options->SetOutputPath(output_path);
}

void BatchReactor::Solve() {
  if (!verbose_) {
    std::cout.setstate(std::ios_base::failbit);
  }

  CeckStatusOfGasMixture();
  // Solve the ODE system: NonIsothermal, Constant Volume
  if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_CONSTANTV) {
    OpenSMOKE::BatchReactor_NonIsothermal_ConstantVolume batch(
        *thermodynamicsMapXML, *kineticsMapXML, *ode_parameters, *batch_options,
        *onTheFlyROPA, *onTheFlyCEMA, *on_the_fly_post_processing, *idt, *polimi_soot,
        volume_, T, P_Pa, omega, global_thermal_exchange_coefficient_, exchange_area_,
        T_environment_);

    batch.Solve(tStart_, tEnd_);
  }

  // Solve the ODE system: NonIsothermal, Volume assigned according to a specified law
  if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_USERDEFINEDVOLUME) {
    std::cout << "No USER DEFINED VOLUME" << std::endl;
    exit(-1);
  }

  // Solve the ODE system: Isothermal, Constant Volume
  if (type_ == OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTV) {
    OpenSMOKE::BatchReactor_Isothermal_ConstantVolume batch(
        *thermodynamicsMapXML, *kineticsMapXML, *ode_parameters, *batch_options,
        *onTheFlyROPA, *onTheFlyCEMA, *on_the_fly_post_processing, *idt, *polimi_soot,
        volume_, T, P_Pa, omega);
    // batch.SaveResults(true);
    batch.Solve(tStart_, tEnd_);

    // Getter
    // time_vector_ = batch.time();
    // temperature_vector_ = batch.final_temperature();
    // pressure_vector_ = batch.final_pressure();
    // mole_fractions_ = batch.mole_fractions();
    // mass_fractions_ = batch.mass_fractions();
  }

  // Solve the ODE system: NonIsothermal, Constant Pressure
  if (type_ == OpenSMOKE::BATCH_REACTOR_NONISOTHERMAL_CONSTANTP) {
    OpenSMOKE::BatchReactor_NonIsothermal_ConstantPressure batch(
        *thermodynamicsMapXML, *kineticsMapXML, *ode_parameters, *batch_options,
        *onTheFlyROPA, *onTheFlyCEMA, *on_the_fly_post_processing, *idt, *polimi_soot,
        volume_, T, P_Pa, omega, global_thermal_exchange_coefficient_, exchange_area_,
        T_environment_);
    batch.Solve(tStart_, tEnd_);
  }

  // Solve the ODE system: Isothermal, Constant Pressure
  if (type_ == OpenSMOKE::BATCH_REACTOR_ISOTHERMAL_CONSTANTP) {
    OpenSMOKE::BatchReactor_Isothermal_ConstantPressure batch(
        *thermodynamicsMapXML, *kineticsMapXML, *ode_parameters, *batch_options,
        *onTheFlyROPA, *onTheFlyCEMA, *on_the_fly_post_processing, *idt, *polimi_soot,
        volume_, T, P_Pa, omega);
    batch.Solve(tStart_, tEnd_);
  }

  OpenSMOKE::OpenSMOKE_logo("OpenSMOKEpp_BatchReactor",
                            "Alberto Cuoci (alberto.cuoci@polimi.it)");
}
