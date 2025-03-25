add_library(OpenSMOKEMaps SHARED OpenSMOKEMaps.cpp)

target_link_libraries(
  OpenSMOKEMaps
  PUBLIC # Eigen
         Eigen3::Eigen
         # Boost
         Boost::date_time
         Boost::filesystem
         Boost::program_options
         Boost::system
         Boost::regex
         Boost::timer
         Boost::chrono
         # pybind
         pybind11::pybind11
         pybind11::python_link_helper
         pybind11::lto
         # Linear Algebra
         ${LINEAR_ALGEBRA_LIB}
)

target_include_directories(
  OpenSMOKEMaps
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}
  PRIVATE ${OPENSMOKE_INCLUDE_DIR}
)

# --------------------------------------------------------------------
# Installing options
# --------------------------------------------------------------------
install(TARGETS OpenSMOKEMaps DESTINATION ${CMAKE_INSTALL_PREFIX})
