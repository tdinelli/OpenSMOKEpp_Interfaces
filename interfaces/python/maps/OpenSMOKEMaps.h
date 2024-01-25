#pragma once

// OpenSMOKE include
#include <OpenSMOKEpp>
#include <maps/Maps_CHEMKIN>

// Standard include
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <string>

namespace py = pybind11;

class OpenSMOKEMaps {
 public:
  OpenSMOKEMaps(const std::string& kinetic_folder, const bool& verbose);

  ~OpenSMOKEMaps();

  void ReadMechanism();

  const OpenSMOKE::ThermodynamicsMap_CHEMKIN* const thermodynamicsMapXML() {
    return thermodynamicsMapXML_;
  }

  const OpenSMOKE::KineticsMap_CHEMKIN* const kineticsMapXML() { return kineticsMapXML_; }

  const static void OpenSMOKEMaps_wrapper(py::module_&);

 private:
  OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML_;
  OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML_;

  boost::filesystem::path kinetics_;
  boost::filesystem::path reaction_names_;
};
