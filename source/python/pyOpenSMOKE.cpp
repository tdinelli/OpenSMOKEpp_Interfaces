// clang-format off
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "maps/OpenSMOKEMaps.h"
#include "ideal_reactors/batch/BatchReactor.h"
// clang-format on

namespace py = pybind11;
constexpr auto call_guard = py::call_guard<py::gil_scoped_release>();

PYBIND11_MODULE(pyOpenSMOKE, m) {
  m.doc() = R"pbdoc(
        pyOpenSMOKE Documentation
        -------------------------

        .. currentmodule:: python_example
        .. autosummary::
            :toctree: _generate

            ThermodynamicsMap_CHEMKIN
            KineticsMap_CHEMKIN
    )pbdoc";

  // Class handler for the maps
  OpenSMOKEMaps::OpenSMOKEMaps_wrapper(m);

  // Thermodynamics Map (Gas Phase)
  // This interface is needed because in this way the module knows the return type of the
  // getter functions that are implemented in OpenSMOKEMaps_wrapper and is able to convert
  // the Object into a python object
  py::class_<OpenSMOKE::ThermodynamicsMap_CHEMKIN>(m, "ThermodynamicsMap_CHEMKIN")
      .def(py::init<const OpenSMOKE::ThermodynamicsMap_CHEMKIN&>(), call_guard,
           "Copy constructor of the original OpenSMOKE class", py::arg("rhs"))
      .def("NumberOfSpecies", &OpenSMOKE::ThermodynamicsMap_CHEMKIN::NumberOfSpecies,
           call_guard,
           "Function that returns the total number of species inside the mechanism")
      .def("NamesOfSpecies", &OpenSMOKE::ThermodynamicsMap_CHEMKIN::NamesOfSpecies,
           call_guard,
           "Function that returns the names of the species within the kinetic mechanism")
      .def("MWs", &OpenSMOKE::ThermodynamicsMap_CHEMKIN::MWs, call_guard,
           "Function that returns a vector containing the molecular weights of the "
           "species [kg/kmol]")
      .def("MW", &OpenSMOKE::ThermodynamicsMap_CHEMKIN::MW, call_guard,
           "Function that returns the molecular weight of the i-th species [kg/kmol]",
           py::arg("i"))
      .def("IndexOfSpecies", &OpenSMOKE::ThermodynamicsMap_CHEMKIN::IndexOfSpecies,
           call_guard,
           "Function that returns the index of the species (1-based) given the name, if the "
           "species is not present returns an error",
           py::arg("name"))
      .def("IndexOfElement", &OpenSMOKE::ThermodynamicsMap_CHEMKIN::IndexOfElement,
           call_guard,
           "Function that returns the index of the element given the name, if the "
           "element is not present returns an error",
           py::arg("name"))
      .def("atomic_composition",
           &OpenSMOKE::ThermodynamicsMap_CHEMKIN::atomic_composition, call_guard,
           "Returns the matrix describing the elemental composition of every species. NS "
           "x NE where NS is number of species and NE is number of elements")
      .def("elements", &OpenSMOKE::ThermodynamicsMap_CHEMKIN::elements, call_guard,
           "Returns the list of names of elements");

  // Kinetic Map (Gas Phase)
  // This interface is needed because in this way the module knows the return type of the
  // getter functions that are implemented in OpenSMOKEMaps_wrapper and is able to convert
  // the Object into a python object
  py::class_<OpenSMOKE::KineticsMap_CHEMKIN>(m, "KineticsMap_CHEMKIN")
      .def(py::init<const OpenSMOKE::KineticsMap_CHEMKIN&>(), call_guard,
           R"pbdoc(Copy constructor of the original OpenSMOKE class)pbdoc",
           py::arg("rhs"))
      .def("SetVerboseOutput", &OpenSMOKE::KineticsMap_CHEMKIN::SetVerboseOutput,
           call_guard, "Sets the verbose output", py::arg("verbose"))
      .def("SetTemperature", &OpenSMOKE::KineticsMap_CHEMKIN::SetTemperature, call_guard,
           "Set the temperature at which the properties have to be evaluated, T in [K]",
           py::arg("T"))
      .def("SetPressure", &OpenSMOKE::KineticsMap_CHEMKIN::SetPressure, call_guard,
           "Set the pressure at which the properties have to be evaluated, P in [Pa]",
           py::arg("P"))
      .def("NamesOfSpecies", &OpenSMOKE::KineticsMap_CHEMKIN::NamesOfSpecies, call_guard,
           "Return the names of all the species")
      .def("FormationRates", &OpenSMOKE::KineticsMap_CHEMKIN::FormationRates, call_guard,
           "Function that computes the formation rates for all the species inside the "
           "kinetic mechanism",
           py::arg("R"))
      .def("HeatRelease", &OpenSMOKE::KineticsMap_CHEMKIN::HeatRelease, call_guard,
           "Function that computes the heat release on the basis of the formation rates "
           "returning a value in [J/m3/s] given as a param the formation rates of all "
           "the species in [kmol/m3/s]",
           py::arg("R"))
      .def("ProductionAndDestructionRates",
           &OpenSMOKE::KineticsMap_CHEMKIN::ProductionAndDestructionRates, call_guard,
           "Calculates the production and the destruction rates for all the species in "
           "the kinetic mechanism [kmol/m3/s], the contributions are calculated on the "
           "basis of net reaction rates",
           py::arg("P"), py::arg("D"))
      .def("ProductionAndDestructionRatesGross",
           &OpenSMOKE::KineticsMap_CHEMKIN::ProductionAndDestructionRatesGross,
           call_guard,
           "Calculates the production and the destruction rates for all the species in "
           "the kinetic mechanism [kmol/m3/s], the contributions are calculated on the "
           "basis of separate forward and backward reaction rates",
           py::arg("P"), py::arg("D"))
      .def("GiveMeReactionRates",
           py::overload_cast<>(&OpenSMOKE::KineticsMap_CHEMKIN::GiveMeReactionRates),
           call_guard, "Returns the net reaction rates in [kmol/m3/s]")
      .def("GiveMeReactionRates",
           py::overload_cast<double*>(
               &OpenSMOKE::KineticsMap_CHEMKIN::GiveMeReactionRates),
           call_guard, "Returns the net reaction rates in [kmol/m3/s]", py::arg("r"))
      .def("GetForwardReactionRates",
           &OpenSMOKE::KineticsMap_CHEMKIN::GetForwardReactionRates, call_guard,
           "Returns the forward reaction rates for all the reactions in the kinetic "
           "scheme [kmol/m3/s]",
           py::arg("r"))
      .def("GetBackwardReactionRates",
           &OpenSMOKE::KineticsMap_CHEMKIN::GetBackwardReactionRates, call_guard,
           "Returns the backward reaction rates for all the reactions in the kinetic "
           "scheme [kmol/m3/s], if a reaction is irreversible, it returns 0",
           py::arg("r"))
      .def("ReactionRates",
           py::overload_cast<const double*>(
               &OpenSMOKE::KineticsMap_CHEMKIN::ReactionRates),
           call_guard,
           "Calculates the reaction rates for all the reactions in the kinetic scheme, c "
           "input param concentration of species in [kmol/m3]",
           py::arg("c"))
      .def("ReactionRates",
           py::overload_cast<const double*, const double>(
               &OpenSMOKE::KineticsMap_CHEMKIN::ReactionRates),
           call_guard,
           "Calculates the reaction rates for all the reactions in the kinetic scheme, c "
           "input param concentration of species in [kmol/m3], cTot total concentration "
           "in [kmol/m3]",
           py::arg("c"), py::arg("cTot"))
      .def("KineticConstants",
           py::overload_cast<>(&OpenSMOKE::KineticsMap_CHEMKIN::KineticConstants),
           call_guard, "Calculates the kinetic constants")
      .def("A", &OpenSMOKE::KineticsMap_CHEMKIN::A, call_guard,
           "Return the frequency factor of a single reaction [kmol, m, s], j index of "
           "reaction 0-based",
           py::arg("j"))
      .def("Beta", &OpenSMOKE::KineticsMap_CHEMKIN::Beta, call_guard,
           "Return the temperature exponent a single reaction, j index of reaction "
           "0-based",
           py::arg("j"))
      .def("E_over_R", &OpenSMOKE::KineticsMap_CHEMKIN::E_over_R, call_guard,
           "Return the activation temperature of a single reaction [K], j index of "
           "reaction 0-based",
           py::arg("j"))
      .def("ThirdBody", &OpenSMOKE::KineticsMap_CHEMKIN::ThirdBody, call_guard,
           "Return the thirdbody efficiency of a species in a given reaction, if a "
           "species does not exist in the list it returns 1. j index of reaction "
           "(0-based), k index of species (0-based)",
           py::arg("j"), py::arg("k"))
      .def("Set_A", &OpenSMOKE::KineticsMap_CHEMKIN::Set_A, call_guard,
           "Sets the frequency factor of a single reaction [kmol, m, s]. j index of the "
           "reaction (0-based), A new value for the frequency factor",
           py::arg("j"), py::arg("A"))
      .def("Set_Beta", &OpenSMOKE::KineticsMap_CHEMKIN::Set_Beta, call_guard,
           "Sets the temperature exponent of a single reaction [K]. j index of the "
           "reaction (0-based), Beta new value for the frequency factor",
           py::arg("j"), py::arg("Beta"))
      .def("Set_EoverR", &OpenSMOKE::KineticsMap_CHEMKIN::Set_E_over_R, call_guard,
           "Sets the activation temperature of a single reaction [K]. j index of the "
           "reaction (0-based), EoverR new value for the frequency factor",
           py::arg("j"), py::arg("E_over_R"))
      .def("Set_ThirdBody", &OpenSMOKE::KineticsMap_CHEMKIN::Set_ThirdBody, call_guard,
           "Sets the third body efficiency of a species in a given reaction. j index of "
           "the reaction (0-based), k index of species (0-based), efficiency new "
           "efficiency coefficient",
           py::arg("j"), py::arg("k"), py::arg("efficiency"))
      .def("kArrheniusModified", &OpenSMOKE::KineticsMap_CHEMKIN::KArrheniusModified,
           call_guard,
           "Returns the modified Arrhenius kinetic constants in [kmol, m, s]");

  // Transport map (Gas Phase) TODO

  // Batch Reactor
  BatchReactor::BatchReactor_wrapper(m);
}
