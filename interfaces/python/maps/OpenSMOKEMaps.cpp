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
