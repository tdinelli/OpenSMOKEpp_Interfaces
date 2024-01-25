#include "OpenSMOKEMaps.h"

#include <boost/filesystem/file_status.hpp>

OpenSMOKEMaps::OpenSMOKEMaps(const std::string& kinetic_folder, const bool& verbose) {
  // TODO: check files and directory existance
  boost::filesystem::path kinetic_folder_ = kinetic_folder;
  kinetics_ = kinetic_folder_ / "kinetics.xml";
  reaction_names_ = kinetic_folder_ / "reaction_names.xml";
};

OpenSMOKEMaps::~OpenSMOKEMaps(){};

void OpenSMOKEMaps::ReadMechanism() {
  // TODO: Implement verbosity
  // Read thermodynamics and kinetics maps
  {
    boost::property_tree::ptree ptree;
    boost::property_tree::read_xml(kinetics_.string(), ptree);

    double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
    thermodynamicsMapXML_ = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(ptree);
    kineticsMapXML_ = new OpenSMOKE::KineticsMap_CHEMKIN(*thermodynamicsMapXML_, ptree);
    double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

    std::cout << "Time to read XML file: " << tEnd - tStart << std::endl;
  }
}

const void OpenSMOKEMaps::OpenSMOKEMaps_wrapper(py::module_& m) {
  constexpr auto call_guard = py::call_guard<py::gil_scoped_release>();
  py::class_<OpenSMOKEMaps>(m, "OpenSMOKEMaps")
      .def(py::init<const std::string&, const bool&>(), call_guard, "Class constructor",
           py::arg("kinetic_folder"), py::arg("verbose"))
      .def("ReadMechanism", &OpenSMOKEMaps::ReadMechanism, call_guard,
           "Function that performs the reading of the mechanism")
      .def("ThermodynamicsMap", &OpenSMOKEMaps::thermodynamicsMapXML, call_guard,
           "Getter function to the smart pointer of the thermodynamic map class")
      .def("KineticsMap", &OpenSMOKEMaps::kineticsMapXML, call_guard,
           "Getter function to the smart pointer of the kinetics map class");
}
