# ----------------------------------------------------------------------------
# Tridiagonal Block matrix decomposition
# ----------------------------------------------------------------------------
set(BLOCKTDSGE_SOURCE_DIR ${PROJECT_SOURCE_DIR}/numerics/blocktdsge)

include(FortranCInterface)
fortrancinterface_verify(CXX)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -m64")

set(BLOCKTDSGE_SOURCE_DIR ${OPENSMOKEPPSOLVERS_INCLUDE_DIR}/numerics/blocktdsge)
add_library(
  blocktdsge
  STATIC
  ${BLOCKTDSGE_SOURCE_DIR}/dgeblttrf.f
  ${BLOCKTDSGE_SOURCE_DIR}/dgeblttrs.f
  ${BLOCKTDSGE_SOURCE_DIR}/aux.f
)
# ----------------------------------------------------------------------------
# Compile definitions to handle additional libraries
# ----------------------------------------------------------------------------
apply_opensmoke_options(blocktdsge DISABLE_ALL)
