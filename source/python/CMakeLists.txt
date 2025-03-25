add_subdirectory(maps)
add_subdirectory(ideal_reactors)

pybind11_add_module(${PROJECT_NAME} ${CMAKE_CURRENT_LIST_DIR}/pyOpenSMOKE.cpp)

target_link_libraries(
  ${PROJECT_NAME}
  PRIVATE Eigen3::Eigen
          Boost::date_time
          Boost::filesystem
          Boost::program_options
          Boost::system
          Boost::regex
          Boost::timer
          Boost::chrono
          # Tridiagonal block solver
          blocktdsge
          # Linear Algebra
          ${LINEAR_ALGEBRA_LIB}
          # Internal modules
          OpenSMOKEMaps
          BatchReactor
)
target_include_directories(
  ${PROJECT_NAME}
  PRIVATE ${CMAKE_CURRENT_LIST_DIR}
          # OpenSMOKEpp
          ${OPENSMOKE_INCLUDE_DIR}
          # Linear algebra acceleration
          ${LINEAR_ALGEBRA_INCLUDE}
)

apply_opensmoke_options(${PROJECT_NAME})

# --------------------------------------------------------------------
# Installing options
# --------------------------------------------------------------------
install(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX})
