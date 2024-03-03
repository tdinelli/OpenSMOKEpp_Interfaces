#--------------------------------------------------------------------
# Preliminary Settings
#--------------------------------------------------------------------

cmake_minimum_required(VERSION 3.20)

project(pyOpenSMOKE
  VERSION 0.1.0
  DESCRIPTION
  "High level interfaces to the OpenSMOKEpp framework"
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

# --------------------------------------------------------------------
# Compulsory (external) libraries
# --------------------------------------------------------------------

# --- Eigen (Compulsory)
find_package(Eigen3 REQUIRED)

# --- Boost (Compulsory)
# --- TODO Clean boost lib I don't think we need all of them
find_package(Boost REQUIRED
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

# --------------------------------------------------------------------
# Subdirectories
# --------------------------------------------------------------------

add_subdirectory(source)

# --------------------------------------------------------------------
# Documentation
# --------------------------------------------------------------------

# TODO

# --------------------------------------------------------------------
# Testing
# --------------------------------------------------------------------

# TODO
