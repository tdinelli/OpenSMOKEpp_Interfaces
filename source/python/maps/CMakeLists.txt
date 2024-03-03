add_library(OpenSMOKEMaps SHARED "")

target_sources(OpenSMOKEMaps
  PRIVATE
  OpenSMOKEMaps.cpp
  PUBLIC
  OpenSMOKEMaps.h
)

target_link_libraries(OpenSMOKEMaps
  PUBLIC
  Eigen3::Eigen
  Boost::date_time
  Boost::filesystem
  Boost::program_options
  Boost::system
  Boost::regex
  Boost::timer
  Boost::chrono
  Boost::serialization
  pybind11::pybind11
  pybind11::python_link_helper
  pybind11::thin_lto
  pybind11::lto
)

target_include_directories(OpenSMOKEMaps
  PUBLIC
  ${CMAKE_CURRENT_LIST_DIR}
  PRIVATE
  ${OPENSMOKE_INCLUDE_DIR}
)

target_compile_options(OpenSMOKEMaps
  PRIVATE
  -Wno-unused-value
  -Wno-mismatched-new-delete
  -Wno-switch
)

# --------------------------------------------------------------------
# Installing options
# --------------------------------------------------------------------

install(TARGETS
  OpenSMOKEMaps
  DESTINATION
  ${CMAKE_INSTALL_PREFIX}
)
