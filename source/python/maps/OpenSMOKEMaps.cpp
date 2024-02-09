#include "OpenSMOKEMaps.h"

OpenSMOKEMaps::OpenSMOKEMaps(const std::string& kinetic_folder, const bool& transport,
                             const bool& verbose) {
  boost::filesystem::path kinetic_folder_ = kinetic_folder;

  kinetics_ = kinetic_folder_ / "kinetics.xml";
  if (!boost::filesystem::exists(kinetics_)) {
    OpenSMOKE::FatalErrorMessage(
        "The folder of the kinetic mechanism does not contains kinetics.xml");
  }

  reaction_names_ = kinetic_folder_ / "reaction_names.xml";
  if (!boost::filesystem::exists(reaction_names_)) {
    OpenSMOKE::FatalErrorMessage(
        "The folder of the kinetic mechanism does not contains reaction_names.xml");
  }
  transport_ = transport;
  verbose_ = verbose;
};

OpenSMOKEMaps::~OpenSMOKEMaps(){};

void OpenSMOKEMaps::ReadMechanism() {
  // Read thermodynamics and kinetics maps
  {
    boost::property_tree::ptree ptree;
    boost::property_tree::read_xml(kinetics_.string(), ptree);

    double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
    thermodynamicsMapXML_ = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(ptree, verbose_);
    kineticsMapXML_ =
        new OpenSMOKE::KineticsMap_CHEMKIN(*thermodynamicsMapXML_, ptree, verbose_);

    if (transport_) {
      transportMapXML_ = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(ptree);
    }
    double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

    if (verbose_) {
      std::cout << "Time to read XML file: " << tEnd - tStart << std::endl;
    }
  }
}

const void OpenSMOKEMaps::OpenSMOKEMaps_wrapper(py::module_& m) {
  // TODO: Improve python side docs.
  constexpr auto call_guard = py::call_guard<py::gil_scoped_release>();
  py::class_<OpenSMOKEMaps>(m, "OpenSMOKEMaps")
      .def(py::init<const std::string&, const bool&, const bool&>(), call_guard,
           "Class constructor", py::arg("kinetic_folder"), py::arg("transport"),
           py::arg("verbose"))
      .def("ReadMechanism", &OpenSMOKEMaps::ReadMechanism, call_guard,
           "Function that performs the reading of the mechanism")
      .def("KineticsMap", &OpenSMOKEMaps::kineticsMapXML, call_guard,
           "Getter function to the smart pointer of the kinetics map class")
      .def("ThermodynamicsMap", &OpenSMOKEMaps::thermodynamicsMapXML, call_guard,
           "Getter function to the smart pointer of the thermodynamic map class")
      .def("TransportMap", &OpenSMOKEMaps::transportMapXML, call_guard,
           "Getter function to the smart pointer of the transport map class");
}
