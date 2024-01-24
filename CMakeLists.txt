#--------------------------------------------------------------------
# Preliminary Settings
#--------------------------------------------------------------------

cmake_minimum_required(VERSION 3.20)

project(
    pyOpenSMOKE
    VERSION 0.1.0
    DESCRIPTION
    "Convenient python interface to OpenSMOKEpp framework"
    LANGUAGES CXX
)

# --- Set module path in order to use custom CMake modules
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_COLOR_MAKEFILE ON)

# --------------------------------------------------------------------
# Compilation options
# --------------------------------------------------------------------

#--------------------------------------------------------------------
# Compulsory (external) libraries
#--------------------------------------------------------------------

# --- Eigen (Compulsory)
find_package(Eigen3 REQUIRED)
# message(STATUS "Eigen libraries: ${EIGEN3_INCLUDE_DIR}")

# --- Boost (Compulsory)
find_package(Boost
    REQUIRED
    COMPONENTS
        date_time
        filesystem
        program_options
        system
        regex
        timer
        chrono
        serialization
)

# --- pybind11 (Compulsory)
find_package(pybind11 REQUIRED)

# --- OpenSMOKE++ (Compulsory)
find_package(OpenSMOKEpp REQUIRED)

# --- OpenSMOKE++ Solvers (Compulsory)
find_package(OpenSMOKEppSolvers REQUIRED)

# --- OpenBLAS library (Compulsory)
set (BLA_VENDOR OpenBLAS)
set (BLA_STATIC ON)
find_package(BLAS REQUIRED)
# list(APPEND LINEAR_ALGEBRA_INCLUDE  $ENV{BLAS_ROOT}/include)
# list(APPEND LINEAR_ALGEBRA_LIB ${BLAS_LIBRARIES})

# --------------------------------------------------------------------
# Sub Modules auto-update (Maybe in the future we will use submodules)
# --------------------------------------------------------------------

# TODO
# include(${CMAKE_MODULE_PATH}/UpdateSubmodules.cmake)

# --------------------------------------------------------------------
# Subdirectories
# --------------------------------------------------------------------

add_subdirectory(interfaces)

# --------------------------------------------------------------------
# Installing options
# --------------------------------------------------------------------

# TODO
# install(TARGETS ${PROJECT_NAME} LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX})

# --------------------------------------------------------------------
# Testing
# --------------------------------------------------------------------

# TODO
