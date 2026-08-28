# Installation

## Supported Platforms

- Linux >=RHEL6
- OS X >=10.5


## Build requirements

- cmake >=3.11 (https://cmake.org)
- cfitsio >=3.430 (https://heasarc.gsfc.nasa.gov/fitsio/fitsio.html)
- gcc >=4.4.7 (https://gcc.gnu.org)
- openmp (http://www.openmp.org/)
- pkg-config (https://www.freedesktop.org/wiki/Software/pkg-config)

## Testing requirements

- python >=3.11 (https://www.python.org)
- pytest
- crds
- ci-watson

GCC may be supplemented by Clang under the following conditions:

- Clang must be compiled with OpenMP support
- A Fortran compiler must be available in `$PATH`


## Prerequisites


### CFITSIO

- Must be compiled as a shared library.

- If `libcfitsio` has been installed to a non-standard path such as `$HOME/programs/cfitsio`, you will need to adjust
  `PKG_CONFIG_PATH` so that `pkg-config` is able to find it:

    ```bash
    export PKG_CONFIG_PATH=$HOME/programs/cfitsio/lib/pkgconfig:$PKG_CONFIG_PATH
    ```

- If you prefer not to use `pkg-config` set `WITH_CFITSIO` to the top-level path where cfitsio is installed.

    ```
    mkdir _build
    cd _build
    cmake .. -DWITH_CFITSIO=$HOME/programs/cfitsio
    ```

### OpenMP

- If you do not have OpenMP or want to disable OpenMP support, set the `ENABLE_OPENMP` option to `OFF`.

    ```
    cmake .. -DENABLE_OPENMP=OFF
    ```

### C unit tests

- If you want to be able to run unit tests written in C.

    ```
    cmake .. -DTESTS=ON
    ```

## Build on Linux

### Debian / Ubuntu

```
apt install cmake libcfitsio-dev gcc gfortran pkg-config
```

### Fedora

```
dnf install cmake cfitsio-devel gcc gcc-gfortran pkgconf-pkg-config
```

### Conda / Mamba

```
conda create -n hstcal -c conda-forge cmake compilers cfitsio pkgconfig python
conda activate hstcal
export LDFLAGS="-Wl,-rpath=$CONDA_PREFIX/lib"
```

### Install

1. Configure

    ```
    mkdir _build
    cd _build
    ```

   If you use Conda/Mamba, you can install it to the env directly; however, note that this would overwrite any previous 
   `hstcal` installation from `conda-forge` channel:

    ```
    cmake .. -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX
    ```

   Otherwise, you can install it to a path you have write access to, e.g.:

    ```
    cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/hstcal
    ```

2. Build

    ```
    make
    ```

3. Install

    ```
    make install
    ```

4. Add HSTCAL to `$PATH` (skip for Conda/Mamba)

    ```bash
    export PATH=$HOME/hstcal/bin:$PATH
    ```

## Build on MacOS / OS X

The LLVM/Clang suite provided by Apple XCode is not sufficient to compile HSTCAL. Please install GCC (C and Fortran
compilers) either from source, or using a package management system such as Homebrew, MacPorts, Fink, or Conda.


### MacPorts

```
port install cmake cfitsio gcc13 pkgconfig +openmp
export CC=gcc-mp-13
export CXX=g++-mp-13
export FC=gfortran-mp-13
export PKG_CONFIG_PATH="/opt/local/lib/pkgconfig:$PKG_CONFIG_PATH"
export LDFLAGS="-Wl,-rpath,/opt/local/lib"
```

### Homebrew

```
brew install cmake cfitsio gcc pkgconfig
export CC=gcc-13
export CXX=gcc-13
export FC=gfortran-13
export PKG_CONFIG_PATH="/opt/homebrew/lib/pkgconfig:$PKG_CONFIG_PATH"
export LDFLAGS="-Wl,-rpath,/opt/homebrew/lib"
```

### Conda / Mamba

```
conda create -n hstcal -c conda-forge cmake compilers cfitsio pkgconfig python
conda activate hstcal
export PKG_CONFIG_PATH="$CONDA_PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH"
export LDFLAGS="-Wl,-rpath,$CONDA_PREFIX/lib"
```

### Install

1. Configure

    ```
    mkdir _build
    cd _build
    ```

   If you use Conda/Mamba, you can install it to the env directly; however, note that this would overwrite any previous
   `hstcal` installation from `conda-forge` channel:

    ```
    cmake .. -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX
    ```

   Otherwise, you can install it to a path you have write access to, e.g.:

    ```
    cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/hstcal
    ```

2. Build

    ```
    make
    ```

3. Install

    ```
    make install
    ```

4. Add HSTCAL to `$PATH` (skip for Conda/Mamba)

    ```bash
    export PATH=$HOME/hstcal/bin:$PATH
    ```

## Debugging

To enable support for debugging symbols use one of the following defines:

```
cmake .. -DCMAKE_BUILD_TYPE=[RelWithDebInfo|Debug]
```

To enable memory leak and heap overflow detection:

```
cmake .. -DENABLE_ASAN=ON [-DENABLE_ASAN_RECOVER=ON]
```

When ASAN (aka AddressAnalyzer) encounters a bug it halts execution and dumps information about the type of error, and
where it occurred. When `ENABLE_ASAN_RECOVER` is enabled, and the `ASAN_OPTIONS` environment variable contains 
`halt_on_error=0`, ASAN will continue to dump information as the program runs. This is mode is incredibly noisy, so it
should only ever be used to test code changes in development.

To show backtraces after calls to `REPORT_ERROR_STATE()` (aka `WhichError`):

```
cmake .. -DENABLE_BACKTRACE=ON
```

Example:

```
# ...
ERROR:    Couldn't process CCD data
ERROR:    CALACS processing NOT completed for file.fits
ERROR:    Invalid temporary file. (status = 1021)
DEBUG:    Location: /path/to/pkg/acs/src/mainacs.c:238
DEBUG:    Backtrace:
/path/to/lib/libhstcalib.so(hstcal_error_state_populate+0xa3) [0x7ddb5c635987]
/path/to/lib/libhstcalib.so(hstcal_error_state_show+0x39) [0x7ddb5c635e8d]
/path/to/pkg/acs/calacs.e(main+0x990) [0x5760d5ad3c4b]
/usr/lib/libc.so.6(+0x27488) [0x7ddb5c035488]
/usr/lib/libc.so.6(__libc_start_main+0x8c) [0x7ddb5c03554c]
/path/to/pkg/acs/calacs.e(_start+0x25) [0x5760d5ad31c5]
```
