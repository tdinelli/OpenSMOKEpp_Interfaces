// pyBIND11 Include
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "ideal_reactors/batch/OpenSMOKE_BatchReactor.h"
#include "OpenSMOKEMaps.h"

namespace py = pybind11;
constexpr auto byref = py::return_value_policy::reference_internal;

PYBIND11_MODULE(OpenSMOKEppInterface, m)
{
	m.doc() = "Python Interface to the Kinetic, Thermodynamics and Transport Maps of OpenSMOKEpp (GAS phase only)";

	// OpenSMOKE Map just a convenient function to properly call
	// a compiled mechanism
	py::class_<OpenSMOKEMaps>(m, "OpenSMOKEMaps")
	.def(py::init<std::string, bool>())
	.def_readwrite("thermodynamicsMapXML", &OpenSMOKEMaps::thermodynamicsMapXML)
	.def_readwrite("kineticsMapXML", &OpenSMOKEMaps::kineticsMapXML)
	;

	// Thermodynamics Map
	py::class_<OpenSMOKE::ThermodynamicsMap_CHEMKIN>(m, "ThermodynamicsMap_CHEMKIN")
	.def("NumberOfSpecies", &OpenSMOKE::ThermodynamicsMap_CHEMKIN::NumberOfSpecies, py::call_guard<py::gil_scoped_release>())
	;
	
	// Kinetic Map
	py::class_<OpenSMOKE::KineticsMap_CHEMKIN>(m, "KineticsMap_CHEMKIN")
	.def("SetTemperature", &OpenSMOKE::KineticsMap_CHEMKIN::SetTemperature, py::call_guard<py::gil_scoped_release>())
	.def("SetPressure", &OpenSMOKE::KineticsMap_CHEMKIN::SetPressure, py::call_guard<py::gil_scoped_release>())
	.def("KineticConstants", py::overload_cast<>(&OpenSMOKE::KineticsMap_CHEMKIN::KineticConstants), py::call_guard<py::gil_scoped_release>())
	.def("kArrheniusModified", &OpenSMOKE::KineticsMap_CHEMKIN::KArrheniusModified, py::call_guard<py::gil_scoped_release>())
	;

	py::class_<BatchReactor>(m, "BatchReactor")
	.def(py::init<std::string>())
	.def("SetInitialComposition", py::overload_cast<std::string, std::vector<std::string>, std::vector<double>>(&BatchReactor::SetInitialComposition), py::call_guard<py::gil_scoped_release>())
	.def("SetTemperature", &BatchReactor::SetTemperature, py::call_guard<py::gil_scoped_release>())
	.def("SetPressure", &BatchReactor::SetPressure, py::call_guard<py::gil_scoped_release>())
	.def("SetDensity", &BatchReactor::SetDensity, py::call_guard<py::gil_scoped_release>())
	.def("SetEndTime", &BatchReactor::SetEndTime, py::call_guard<py::gil_scoped_release>())
	.def("SetStartTime", &BatchReactor::SetStartTime, py::call_guard<py::gil_scoped_release>())
	.def("SetVolume", &BatchReactor::SetVolume, py::call_guard<py::gil_scoped_release>())
	.def("SetExchangeArea", &BatchReactor::SetExchangeArea, py::call_guard<py::gil_scoped_release>())
	.def("Set_global_thermal_exchange_coefficient", &BatchReactor::Set_global_thermal_exchange_coefficient, py::call_guard<py::gil_scoped_release>())
	.def("SetEnvironmentTemperature", &BatchReactor::SetEnvironmentTemperature, py::call_guard<py::gil_scoped_release>())
	.def("SetType", &BatchReactor::SetType, py::call_guard<py::gil_scoped_release>())
	.def("Solve", &BatchReactor::Solve, py::call_guard<py::gil_scoped_release>())
	.def_property_readonly("tempo", &BatchReactor::GetTempo)
	.def_property_readonly("temperature", &BatchReactor::GetTemperatura)
	.def_property_readonly("pressure", &BatchReactor::GetPressione)
	.def("mole_fractions", &BatchReactor::GetMoles)
	;
	// .def("masses", &BatchReactor::GetMasses)

}
