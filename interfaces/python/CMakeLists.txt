add_subdirectory(maps)
add_subdirectory(ideal_reactors)

pybind11_add_module(${PROJECT_NAME}
    ${CMAKE_CURRENT_LIST_DIR}/pyOpenSMOKE.cpp
)

target_link_libraries(${PROJECT_NAME}
    PRIVATE
        Eigen3::Eigen
        Boost::date_time
        Boost::filesystem
        Boost::program_options
        Boost::system
        Boost::regex
        Boost::timer
        Boost::chrono
        Boost::serialization
        OpenSMOKEMaps
)

target_include_directories(${PROJECT_NAME}
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}
    PRIVATE
        ${OPENSMOKE_INCLUDE_DIR}
)

target_compile_options(${PROJECT_NAME}
    PRIVATE
        -Wno-unused-value
        -Wno-mismatched-new-delete
        -Wno-switch
)
