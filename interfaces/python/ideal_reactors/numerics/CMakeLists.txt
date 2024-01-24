project (blocktdsge
    DESCRIPTION
    "TODO"
    LANGUAGES 
    C CXX Fortran)

# --- Set module path in order to use custom CMake modules
set (PROJECT_SOURCE_DIR ${OPENSMOKEPPSOLVERS_INCLUDE_DIR}/numerics/blocktdsge)

# --- Controlling fortran c++ interoperability
include(FortranCInterface)
FortranCInterface_VERIFY(CXX)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -m64")

add_library(blocktdsge
    ${PROJECT_SOURCE_DIR}/dgeblttrf.f
    ${PROJECT_SOURCE_DIR}/dgeblttrs.f
    ${PROJECT_SOURCE_DIR}/aux.f
)
