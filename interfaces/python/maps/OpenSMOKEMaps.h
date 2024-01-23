#pragma once

// OpenSMOKE include
#include <OpenSMOKEpp>
#include <maps/Maps_CHEMKIN>

// Standard include
#include <string>

class OpenSMOKEMaps {
 public:
  OpenSMOKEMaps(const std::string& kinetic_folder, const bool& verbose);

  ~OpenSMOKEMaps();

  void ReadMechanism();

  const OpenSMOKE::ThermodynamicsMap_CHEMKIN* const thermodynamicsMapXML() {
    return thermodynamicsMapXML_;
  }

  const OpenSMOKE::KineticsMap_CHEMKIN* const kineticsMapXML() { return kineticsMapXML_; }

 private:
  OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML_;
  OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML_;

  boost::filesystem::path kinetics_;
  boost::filesystem::path reaction_names_;
};
