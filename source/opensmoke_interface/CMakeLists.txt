#--------------------------------------------------------------------
# Preliminary Settings
#--------------------------------------------------------------------

cmake_minimum_required (VERSION 3.20)

project (OpenSMOKEppInterface)

# --- Set module path in order to use custom CMake modules
set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

set (CMAKE_CXX_STANDARD 17)
set (PYBIND11_PYTHON_VERSION 3.11)
set (CMAKE_POSITION_INDEPENDENT_CODE ON)
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{EIGEN_LIBRARY_PATH}/cmake)
set (CMAKE_COLOR_MAKEFILE ON)
set (VERSION "0.2.0")
set (INSTALL_BIN_DIR "${CMAKE_INSTALL_PREFIX}/bin" CACHE PATH "Installation directory for executables")
#set (MACOSX_DEPLOYMENT_TARGET 12.6)

if( NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release)
endif()

message (STATUS "Compiling in mode: ${CMAKE_BUILD_TYPE}")

#--------------------------------------------------------------------
# Compulsory libraries
#--------------------------------------------------------------------

# --- Source code for the project
message (STATUS "CurveMatching source: ${PROJECT_SOURCE_DIR}")

# --- Eigen (Compulsory)
find_package (Eigen3 REQUIRED)
message (STATUS "Eigen libraries: ${EIGEN3_INCLUDE_DIR}")

# --- Boost (Compulsory)
find_package (Boost REQUIRED COMPONENTS date_time filesystem program_options system regex timer chrono serialization)
message (STATUS "Boost include directory: ${Boost_INCLUDE_DIRS}")
message (STATUS "Boost libraries: ${Boost_LIBRARIES}")

# --- OpenSMOKE++ (Compulsory)
find_package(OpenSMOKEpp REQUIRED)
message (STATUS "OpenSMOKE++: ${OpenSMOKEpp_INCLUDE_DIR}")

# --- pybind11 (Compulsory)
find_package(pybind11 REQUIRED)
message(STATUS "pyBIND11: ${PYBIND11_INCLUDE_DIRS}")

set(COMPULSORY_LIBRARIES ${Boost_LIBRARIES})

set(COMPULSORY_INCLUDE ${Boost_INCLUDE_DIRS} 
                       ${PYBIND11_INCLUDE_DIRS}
                       ${EIGEN3_INCLUDE_DIR}
                       ${OpenSMOKEpp_INCLUDE_DIR})

add_compile_options(-std=c++17 -fPIC -m64 -w -Wall -Wextra)
add_link_options(-std=c++17 -fPIC -m64 -w -Wall -Wextra)

#--------------------------------------------------------------------
# Include flags
#--------------------------------------------------------------------

include_directories (${PROJECT_SOURCE_DIR}
                     ${COMPULSORY_INCLUDE})

#--------------------------------------------------------------------
# Library flags
#--------------------------------------------------------------------

link_libraries (${COMPULSORY_LIBRARIES})

#--------------------------------------------------------------------
# Targets to compile
#--------------------------------------------------------------------

pybind11_add_module(
	${PROJECT_NAME}
	opensmoke_interface/OpenSMOKEppInterface_py.cpp
)

#target_compile_definitions(${PROJECT_NAME})# PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
target_include_directories(${PROJECT_NAME} PRIVATE ${PYBIND11_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR})
target_link_libraries(${PROJECT_NAME} PRIVATE Eigen3::Eigen)

#--------------------------------------------------------------------
# Installing options
#--------------------------------------------------------------------

#install(TARGETS bspl.sh DESTINATION ${INSTALL_BIN_DIR})

#--------------------------------------------------------------------
# Testing options
#--------------------------------------------------------------------

# enable_testing()
