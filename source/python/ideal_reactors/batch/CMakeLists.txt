add_library(BatchReactor SHARED BatchReactor.cpp)

target_link_libraries(
  BatchReactor
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
         # BLOCKTDSOLVER
         blocktdsge
         # Linear Algebra
         ${LINEAR_ALGEBRA_LIB}
         # pybind
         pybind11::pybind11
         pybind11::python_link_helper
         pybind11::lto
)

target_include_directories(
  BatchReactor PUBLIC ${CMAKE_CURRENT_LIST_DIR} ${OPENSMOKE_INCLUDE_DIR} ${LINEAR_ALGEBRA_INCLUDE}
)

apply_opensmoke_options(BatchReactor)

# --------------------------------------------------------------------
# Installing options
# --------------------------------------------------------------------
install(TARGETS BatchReactor DESTINATION ${CMAKE_INSTALL_PREFIX})
