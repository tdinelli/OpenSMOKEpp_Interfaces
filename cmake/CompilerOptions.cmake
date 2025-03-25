# ----------------------------------------------------------------------------
# Determine System Architecture (32/64 bit)
# ----------------------------------------------------------------------------
set(X64 OFF)
if(CMAKE_SIZEOF_VOID_P EQUAL 8)
  set(X64 ON)
endif()
message(STATUS "Architecture X64: " ${X64})

# ----------------------------------------------------------------------------
# GNU or CLANG compilers
# ----------------------------------------------------------------------------
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU" OR "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")

  # --------------------------------------------------------------------------
  # Set C++ Standard
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF) # Ensures that no GNU extensions are used

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-value")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-mismatched-new-delete")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-switch")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-dangling-else")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-return-type-c-linkage")

  # --------------------------------------------------------------------------
  # Debug options
  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g --coverage -fprofile-arcs -ftest-coverage")

  # --------------------------------------------------------------------------
  # Release options
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fPIC -march=native")
endif()
