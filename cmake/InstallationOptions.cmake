# --------------------------------------------------------------------
# Installation directories
# --------------------------------------------------------------------
include(GNUInstallDirs)

# --------------------------------------------------------------------
# Define installation directories
set(INSTALL_BIN_DIR
    ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}
    CACHE PATH "Installation directory for executables")
set(INSTALL_LIB_DIR
    ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
    CACHE PATH "Installation directory for libraries")

# --------------------------------------------------------------------
# Function to handle automatically executable installation
function(install_executable_with_deps target_name)
  install(
    TARGETS ${target_name}
    RUNTIME DESTINATION ${INSTALL_BIN_DIR}
    LIBRARY DESTINATION ${INSTALL_LIB_DIR}
    ARCHIVE DESTINATION ${INSTALL_LIB_DIR})

  install(
    IMPORTED_RUNTIME_ARTIFACTS
    Boost::date_time
    Boost::filesystem
    Boost::program_options
    Boost::system
    Boost::regex
    Boost::timer
    Boost::chrono
    RUNTIME DESTINATION ${INSTALL_BIN_DIR}
    LIBRARY DESTINATION ${INSTALL_LIB_DIR}
  )

  install(
    CODE "file(GET_RUNTIME_DEPENDENCIES
    EXECUTABLES \$<TARGET_FILE:${target_name}>
    RESOLVED_DEPENDENCIES_VAR _r_deps
    UNRESOLVED_DEPENDENCIES_VAR _u_deps
    PRE_INCLUDE_REGEXES \".*\"
    PRE_EXCLUDE_REGEXES \"^/System/.*\" \"^/usr/lib/.*\"
  )
  foreach(_file \${_r_deps})
    file(INSTALL
      DESTINATION \"${INSTALL_LIB_DIR}\"
      TYPE SHARED_LIBRARY
      FILES \"\${_file}\"
    )
  endforeach()
  list(LENGTH _u_deps _u_deps_length)
  if(_u_deps_length GREATER 0)
    message(STATUS \"Unresolved dependencies: \${_u_deps}\")
  endif()")
endfunction()
