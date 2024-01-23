#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "maps/OpenSMOKEMaps.h"

namespace py = pybind11;

constexpr auto byref = py::return_value_policy::reference_internal;
constexpr auto call_guard = py::call_guard<py::gil_scoped_release>();

PYBIND11_MODULE(pyOpenSMOKE, m) {
  m.doc() = "TODO";

  // Class handler for the maps
  py::class_<OpenSMOKEMaps>(m, "OpenSMOKEMaps")
      .def(py::init<const std::string&, const bool&>(), call_guard, "Class constructor",
           py::arg("kinetic_folder"), py::arg("verbose"))
      .def("ReadMechanism", &OpenSMOKEMaps::ReadMechanism, call_guard,
           "Function that performs the reading of the mechanism")
      .def("ThermodynamicsMap", &OpenSMOKEMaps::thermodynamicsMapXML, call_guard,
           "Getter function to the smart pointer of the thermodynamic map class")
      .def("KineticsMap", &OpenSMOKEMaps::kineticsMapXML, call_guard,
           "Getter function to the smart pointer of the kinetics map class");

  // Thermodynamics Map (Gas Phase)
  // This interface is needed because in this way the module knows the return type and is
  // able to convert the Object into a python object
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
           py::arg("i"));

  // Kinetic Map (Gas Phase)
  // This interface is needed because in this way the module knows the return type and is
  // able to convert the Object into a python object
  py::class_<OpenSMOKE::KineticsMap_CHEMKIN>(m, "KineticsMap_CHEMKIN")
      .def(py::init<const OpenSMOKE::KineticsMap_CHEMKIN&>(), call_guard,
           "Copy constructor of the original OpenSMOKE class", py::arg("rhs"))
      .def("SetTemperature", &OpenSMOKE::KineticsMap_CHEMKIN::SetTemperature, call_guard,
           "Set the temperature at which the properties have to be evaluated, T in [K]",
           py::arg("T"))
      .def("SetPressure", &OpenSMOKE::KineticsMap_CHEMKIN::SetPressure, call_guard,
           "Set the pressure at which the properties have to be evaluated, P in [Pa]",
           py::arg("P"))
      .def("NamesOfSpecies", &OpenSMOKE::KineticsMap_CHEMKIN::NamesOfSpecies, call_guard,
           "Return the namse of the species")
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
      .def("KineticConstants",
           py::overload_cast<>(&OpenSMOKE::KineticsMap_CHEMKIN::KineticConstants),
           call_guard, "Function that calculates the kinetic constants")
      .def("kArrheniusModified", &OpenSMOKE::KineticsMap_CHEMKIN::KArrheniusModified,
           call_guard,
           "Returns the modified Arrhenius kinetic constants in [kmol, m, s]");

    // Transport map (Gas Phase) TODO

    // Batch Reactor
}
