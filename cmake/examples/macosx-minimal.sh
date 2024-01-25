cmake \
    -D CMAKE_CXX_COMPILER:PATH=/opt/homebrew/opt/llvm/bin/clang++ \
    -D CMAKE_EXPORT_COMPILE_COMMANDS:BOOL=TRUE \
    -D Eigen3_DIR:PATH=/Users/tdinelli/NumericalLibraries/eigen/eigen-3.4.0/share/eigen3/cmake \
    -D Boost_ROOT:PATH=/Users/tdinelli/NumericalLibraries/boost/boost-1.83.0-clang-17.0.1 \
    -D OpenSMOKEpp_ROOT:PATH=/Users/tdinelli/Documents/GitHub/OpenSMOKEpp \
    -D OpenSMOKEppSolvers_ROOT:PATH=/Users/tdinelli/Documents/GitHub/OpenSMOKEppSolvers \
    -D BLAS_ROOT:PATH=/Users/tdinelli/NumericalLibraries/openblas/openblas-0.3.24-clang-17.0.3 \
    ..
