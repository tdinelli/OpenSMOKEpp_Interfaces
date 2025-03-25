# ----------------------------------------------------------------------------
# List containing the options activated during the configuration of the CMake project
# ----------------------------------------------------------------------------
set(PROGRAM_OPTIONS)

# ----------------------------------------------------------------------------
# Function to apply OpenSMOKE options to a specific target
function(apply_opensmoke_options TARGET_NAME)
  # Parse function arguments
  set(options DISABLE_ALL) # Option to disable all features
  set(oneValueArgs "") # Arguments that take one value
  set(multiValueArgs OPTIONS) # Arguments that take multiple values
  cmake_parse_arguments(
    ARG
    "${options}"
    "${oneValueArgs}"
    "${multiValueArgs}"
    ${ARGN}
  )

  # Get the list of options to process
  if(NOT ARG_OPTIONS)
    # Default options list if no custom options provided
    set(process_options
        OPENSMOKE_USE_CEQ
        OPENSMOKE_USE_MKL
        OPENSMOKE_USE_OPENBLAS
        OPENSMOKE_USE_SUPERLU_SERIAL
        OPENSMOKE_USE_LINPACK
        OPENSMOKE_USE_SUNDIALS
        OPENSMOKE_USE_UMFPACK
        OPENSMOKE_USE_BZZMATH
        OPENSMOKE_USE_RADAU
        OPENSMOKE_USE_ODEPACK
        OPENSMOKE_USE_DASPK
        OPENSMOKE_USE_MEBDF
    )
  else()
    set(process_options ${ARG_OPTIONS})
  endif()

  foreach(OPTION IN LISTS process_options)
    if(ARG_DISABLE_ALL)
      target_compile_definitions(${TARGET_NAME} PRIVATE ${OPTION}=0)
    else()
      if(${OPTION})
        target_compile_definitions(${TARGET_NAME} PRIVATE ${OPTION}=1)
      else()
        target_compile_definitions(${TARGET_NAME} PRIVATE ${OPTION}=0)
      endif()
    endif()
  endforeach()
endfunction()

# ----------------------------------------------------------------------------
# Options
# ----------------------------------------------------------------------------
option(OPENSMOKE_USE_OPENMP "Activate support for OpenMP parallelization" OFF)
if(OPENSMOKE_USE_OPENMP AND NOT OPENSMOKE_DISABLE_ALL)
  include(${CMAKE_CURRENT_LIST_DIR}/OpenMP.cmake)
endif()

# ----------------------------------------------------------------------------
# SUPPORT FOR THERMODYNAMIC EQUILIBRIUM CALCULATION THROUGH CEQ LIBRARY
# https://tcg.mae.cornell.edu/CEQ/
# ----------------------------------------------------------------------------
option(OPENSMOKE_USE_CEQ
       "Activate CEQ support for OpenSMOKEpp to handle thermodynamic equilibrium calculations" OFF
)

# ----------------------------------------------------------------------------
# Libconfig library from https://hyperrealm.github.io/libconfig/ Needed specifically by the
# CounterFlowFlame and MicrogravityDroplet solvers
# ----------------------------------------------------------------------------
option(OPENSMOKE_USE_CONFIG "Activate support for libconfig" OFF)

# ----------------------------------------------------------------------------
# Linear Algebra Acceleration
# ----------------------------------------------------------------------------
option(OPENSMOKE_USE_MKL "Activate Intel MKL linear algebra accelaration for OpenSMOKEpp" OFF)
option(OPENSMOKE_USE_OPENBLAS "Activate OpenBLAS linear algebra acceleration for OpenSMOKEpp" OFF)
if(OPENSMOKE_USE_MKL
   AND OPENSMOKE_USE_OPENBLAS
   AND NOT OPENSMOKE_DISABLE_ALL
)
  message(
    FATAL_ERROR
      "Solvers can be compiled with just one BLAS/LAPACK distribution at a time! Choose between MKL and OpenBLAS"
  )
endif()

option(OPENSMOKE_USE_SUPERLU_SERIAL "Activate SuperLU-serial support for OpenSMOKEpp" OFF)
option(OPENSMOKE_USE_LINPACK "Activate LINPACK support for OpenSMOKEpp" OFF)

# ----------------------------------------------------------------------------
# Additional ODE and DAE Solvers libraries
# ----------------------------------------------------------------------------
option(OPENSMOKE_USE_SUNDIALS "Activate SUNDIALS support for OpenSMOKEpp" OFF)
option(OPENSMOKE_USE_UMFPACK "Activate UMFPACK support for OpenSMOKEpp" OFF)
option(OPENSMOKE_USE_BZZMATH "Activate BzzMath support for OpenSMOKEpp" OFF)
option(OPENSMOKE_USE_RADAU "Activate RADAU support for OpenSMOKEpp" OFF)
option(OPENSMOKE_USE_ODEPACK "Activate ODEPACK support for OpenSMOKEpp" OFF)
option(OPENSMOKE_USE_DASPK "Activate DASPK support for OpenSMOKEpp" OFF)
option(OPENSMOKE_USE_MEBDF "Activate MEBDF support for OpenSMOKEpp" OFF)
