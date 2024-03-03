# OpenSMOKEpp (High Level) interfaces

Repository containing my high level programming language interfaces (e.g. [python](source/python), Julia (one day)) for OpenSMOKEpp.

## OpenSMOKEpp interfaced classes

- Gas Phase kinetic map.
- Gas Phase Thermodynamic properties map.
- Gas Phase Transport properties map. (**WIP**)
- Batch Reactor.

## EXAMPLES
**Incoming**

## INSTALLATION

**Disclaimer** the package is still under active development for this reason at the moment has been
tested only on my personal PC and in some other controlled environments, it is not guaranteed to
work as is.

### Requirements
- **CMake** (https://cmake.org/).
- **Ninja** (https://ninja-build.org/).
- **OpenSMOKEpp** (https://www.opensmokepp.polimi.it/) distributed under request contact Alberto Cuoci <alberto.cuoci@polimi.it>.
- **OpenSMOKEppSolvers** (https://www.opensmokepp.polimi.it/) distributed under request contact Alberto Cuoci <alberto.cuoci@polimi.it>.
- **Eigen3** (https://eigen.tuxfamily.org/).
- **Boost** (https://www.boost.org/).
- **OpenBLAS** (https://www.openblas.net/).
- Working C/C++ compiler (e.g. [gnu](https://gcc.gnu.org/) or [llvm](https://llvm.org/) compiler).
- Working Fortran compiler ([gfortran](https://gcc.gnu.org/wiki/GFortran)).
- Strongly recommended Anaconda/Miniconda (https://anaconda.org/).

Once all the necessary requirements are satisfied in order to build and install the package follow the following instructions:

1. Clone this repository into your system.
    ```bash
    > git clone https://github.com/tdinelli/OpenSMOKEpp_Interfaces.git --depth=1
    ```
2. Move into the main folder.
    ```bash
    > cd OpenSMOKEpp_Interfaces
    ```
3. Export the following environment variables (with the appropriate paths).
    ```bash
    > export CC="/opt/homebrew/opt/llvm/bin/clang"
    > export CXX="/opt/homebrew/opt/llvm/bin/clang++"
    > export FC="/opt/homebrew/bin/gfortran-13"
    > export Boost_ROOT="$HOME/NumericalLibraries/boost/boost-1.83.0-clang-17.0.1"
    > export BLAS_ROOT="$HOME/NumericalLibraries/openblas/openblas-0.3.24-clang-17.0.3"
    > export OpenSMOKEpp_ROOT="$HOME/Documents/GitHub/OpenSMOKEpp"
    > export OpenSMOKEppSolvers_ROOT="$HOME/Documents/GitHub/OpenSMOKEppSolvers"
    > export Eigen3_ROOT="$HOME/NumericalLibraries/eigen/eigen-3.4.0/share/eigen3/cmake"
    ```
4. Run the following commands. In order to avoid problems with python versioning and to avoid conflicts between existing python packages it is highly recommended to install the package into a (clean) conda environment.
    ```bash
    > python3 -m pip install .
    ```
    If you need to debug the installation run the previous command with "**-v**" at the end. If everything is done correct something like this should appear on the screen.
    ```bash
    ...
    Successfully built pyOpenSMOKE
    Installing collected packages: pyOpenSMOKE
    Successfully installed pyOpenSMOKE-0.0.1
    ```
