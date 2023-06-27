// pyBIND11 Include
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

// Eigen Include
#include <Eigen/Dense>
#include <Eigen/Sparse>

// OpenSMOKEpp Include
#include "OpenSMOKEpp"
#include "maps/Maps_CHEMKIN"

namespace py = pybind11;
constexpr auto byref = py::return_value_policy::reference_internal;

struct KineticMap
{
	OpenSMOKE::ThermodynamicsMap_CHEMKIN*	thermodynamicsMapXML;
	OpenSMOKE::KineticsMap_CHEMKIN*			kineticsMapXML;

	KineticMap(std::string path_kinetics, bool verbose)
	{
		boost::filesystem::path path_kinetics_ = path_kinetics;

		if (!verbose) std::cout.setstate(std::ios_base::failbit);
		
		{
			// Read thermodynamics and kinetics maps
			boost::property_tree::ptree ptree;
    		boost::property_tree::read_xml( (path_kinetics_ / "kinetics.xml").string(), ptree );

			double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			thermodynamicsMapXML = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(ptree);
			kineticsMapXML = new OpenSMOKE::KineticsMap_CHEMKIN(*thermodynamicsMapXML, ptree);
			double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

			std::cout << "Time to read XML file: " << tEnd - tStart << std::endl;
		}
		std::cout.clear();
	}
};

PYBIND11_MODULE(OpenSMOKEppInterface, m)
{
	m.doc() = "Python Interface to the Kinetic, Thermodynamics and Transport Maps of OpenSMOKEpp (GAS phase only)";

	py::class_<KineticMap>(m, "KineticMap")
	.def(py::init<std::string, bool>())
	.def_readwrite("thermodynamicsMapXML", &KineticMap::thermodynamicsMapXML)
	.def_readwrite("kineticsMapXML", &KineticMap::kineticsMapXML)
	;

	py::class_<OpenSMOKE::ThermodynamicsMap_CHEMKIN>(m, "ThermodynamicsMap_CHEMKIN")
	.def("NumberOfSpecies", &OpenSMOKE::ThermodynamicsMap_CHEMKIN::NumberOfSpecies, py::call_guard<py::gil_scoped_release>())
	;
	
	py::class_<OpenSMOKE::KineticsMap_CHEMKIN>(m, "KineticsMap_CHEMKIN")
	.def("SetTemperature", &OpenSMOKE::KineticsMap_CHEMKIN::SetTemperature, py::call_guard<py::gil_scoped_release>())
	.def("SetPressure", &OpenSMOKE::KineticsMap_CHEMKIN::SetPressure, py::call_guard<py::gil_scoped_release>())
	.def("KineticConstants", py::overload_cast<>(&OpenSMOKE::KineticsMap_CHEMKIN::KineticConstants), py::call_guard<py::gil_scoped_release>())
	.def("kArrheniusModified", &OpenSMOKE::KineticsMap_CHEMKIN::KArrheniusModified, py::call_guard<py::gil_scoped_release>())
	;
}
