add_subdirectory(maps)
add_subdirectory(ideal_reactors)

pybind11_add_module(${PROJECT_NAME}
    ${CMAKE_CURRENT_LIST_DIR}/pyOpenSMOKE.cpp
)

target_link_libraries(${PROJECT_NAME}
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
    BLAS::BLAS
    blocktdsge
    OpenSMOKEMaps
    BatchReactor
)

target_include_directories(${PROJECT_NAME}
    PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
    ${OPENSMOKE_INCLUDE_DIR}
    ${BLAS_ROOT}/include
)

target_compile_options(${PROJECT_NAME}
    PRIVATE
    -Wno-unused-value
    -Wno-mismatched-new-delete
    -Wno-switch
    -Wno-dangling-else
    -Wno-return-type-c-linkage
)

target_compile_definitions(${PROJECT_NAME}
    PRIVATE
    OPENSMOKE_USE_OPENBLAS=1
    VERSION_INFO=${EXAMPLE_VERSION_INFO}
)

# --------------------------------------------------------------------
# Installing options
# --------------------------------------------------------------------

install(TARGETS
    ${PROJECT_NAME}
    DESTINATION
    ${CMAKE_INSTALL_PREFIX}
)
