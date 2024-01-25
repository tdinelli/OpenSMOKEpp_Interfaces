add_library(BatchReactor SHARED "")

target_sources(BatchReactor
    PRIVATE
        BatchReactor.cpp
    PUBLIC
        BatchReactor.h
)

target_link_libraries(BatchReactor
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
        pybind11::pybind11
        pybind11::python_link_helper
        pybind11::thin_lto
        pybind11::lto
    PRIVATE
        blocktdsge
)

target_include_directories(BatchReactor
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}
    PRIVATE
        ${OPENSMOKE_INCLUDE_DIR}
        ${BLAS_ROOT}/include
)

target_compile_options(BatchReactor
    PRIVATE
        -Wno-unused-value
        -Wno-mismatched-new-delete
        -Wno-switch
        -Wno-self-assign-field
        -Wno-dangling-else
        -Wno-return-type-c-linkage
)

target_compile_definitions(BatchReactor
    PRIVATE
        OPENSMOKE_USE_OPENBLAS=1
)
