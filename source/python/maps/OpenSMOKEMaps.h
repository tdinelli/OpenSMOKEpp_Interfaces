#pragma once

#include <OpenSMOKEpp>
#include <maps/Maps_CHEMKIN>

#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

/**
 * @class OpenSMOKEMaps
 * @brief Class that handle the reading of the kinetic mechanism, provided in .XML format
 * coherent to the latest version of OpenSMOKEpp. The implementation of this class can
 * seem useless but since the base constructor of the kinetic and thermo map of
 * OpenSMOKEpp accept as input parameter the pointer to the boost ptree XML reader this is
 * why this class exists, on top of that there is a copy constructor that accept a raw
 * pointer so I cannot use smart pointer in this context, which is fine at the moment but
 * this behaviour has to change in the future. Further details on how to use smart
 * pointers in pybind can be found here
 * https://pybind11.readthedocs.io/en/stable/advanced/smart_ptrs.html
 */
class OpenSMOKEMaps {
 public:
  /**
   * @brief Default constructor for the class. Some consistency checks are performed here.
   *
   * @param kinetic_folder string representing the path to the folder containing the
   * reaction_names.xml and kinetics.xml files.
   * @param transport boolean value that activate the generation of the transport
   * properties map
   * @param verbose boolean value to set the logging of relevant informations when parsing
   * the mechanism.
   */
  OpenSMOKEMaps(const std::string& kinetic_folder, const bool& transport, const bool& verbose);

  /**
   * @brief Default destructor of the class.
   */
  ~OpenSMOKEMaps();

  /**
   * @brief Function that performs the actual reading of the kinetic scheme.
   */
  void ReadMechanism();

  /**
   * @brief Getter function pointing to the object representing the thermodynamics map.
   *
   * @return Pointer to the thermodynamic map.
   */
  const OpenSMOKE::ThermodynamicsMap_CHEMKIN* const thermodynamicsMapXML() {
    return thermodynamicsMapXML_;
  }

  /**
   * @brief Getter function pointing to the object representing the kinetics map.
   *
   * @return Pointer to the kinetics map.
   */
  const OpenSMOKE::KineticsMap_CHEMKIN* const kineticsMapXML() {
    return kineticsMapXML_;
  }

  /**
   * @brief Getter function pointing to the object representing the transport map.
   *
   * @return Pointer to the transport map.
   */
  const OpenSMOKE::TransportPropertiesMap_CHEMKIN* const transportMapXML() {
    return transportMapXML_;
  }

  /**
   * @brief Function that handle the wrap to the python objects. This is here in order to
   * speed up the compilation process. For reference see:
   * https://pybind11.readthedocs.io/en/stable/advanced/misc.html#partitioning-code-over-multiple-extension-modules
   */
  const static void OpenSMOKEMaps_wrapper(py::module_&);

 private:
  OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML_;  //!< Raw pointer to the thermodynamic map.

  OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML_;  //!< Raw pointer to the kinetic map.

  OpenSMOKE::TransportPropertiesMap_CHEMKIN* transportMapXML_;  //!< Raw pointer to the transport map.

  boost::filesystem::path kinetics_;  //!< Path to the kinetics.xml file.

  boost::filesystem::path reaction_names_;  //!< Path to the reaction_names.xml file.

  bool verbose_;  //!< Variable to control the verbosity of the reading.

  bool transport_;  //!< Whether to process or not the transport map.
};
