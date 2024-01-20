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

- python >=3.8 (https://www.python.org)

GCC may be supplemented by Clang under the following conditions:

- Clang must be compiled with OpenMP support
- A Fortran compiler must be available in `$PATH`


## Prerequisites


### CFITSIO

- Must be compiled as a shared library.

- If `libcfitsio` has been installed to a non-standard path such as `$HOME/programs/cfitsio`, you will need to adjust `PKG_CONFIG_PATH` so that `pkg-config` is able to find it:

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
export PKG_CONFIG_PATH="$CONDA_PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH"
export LDFLAGS="-Wl,-rpath=$CONDA_PREFIX/lib"
```

1. Configure

    ```
    mkdir _build
    cd _build
    cmake .. -DCMAKE_INSTAL_PREFIX=$HOME/hstcal
    ```

2. Build

    ```
    make
    ```

3. Install

    ```
    make install
    ```

4. Add HSTCAL to `$PATH`

    ```bash
    export PATH=$HOME/hstcal/bin:$PATH
    ```

## Build on MacOS / OS X

The LLVM/Clang suite provided by Apple XCode is not sufficient to compile HSTCAL. Please install GCC (C and Fortran compilers) either from source, or using a package management system such as Homebrew, MacPorts, Fink, or Conda.


### MacPorts

```
port install cmake cfitsio gcc13 pkgconfig +openmp
export PKG_CONFIG_PATH="/opt/local/lib/pkgconfig:$PKG_CONFIG_PATH"
export LDFLAGS="-Wl,-rpath,/opt/local/lib"
```

### Homebrew

```
brew install cmake cfitsio gcc pkgconfig
export PKG_CONFIG_PATH="/opt/homebrew/lib/pkgconfig:$PKG_CONFIG_PATH"
export LDFLAGS="-Wl,-rpath,/opt/homebrew/lib"
```

### Fink

```
fink install cmake libcfitsio10-dev gcc12 pkgconfig
export PKG_CONFIG_PATH="/opt/sw/lib/pkgconfig:$PKG_CONFIG_PATH"
export LDFLAGS="-Wl,-rpath,/opt/sw/lib"
```

### Conda / Mamba

```
conda create -n hstcal -c conda-forge cmake compilers cfitsio pkgconfig python
conda activate hstcal
export PKG_CONFIG_PATH="$CONDA_PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH"
export LDFLAGS="-Wl,-rpath,$CONDA_PREFIX/lib"
```

1. Configure

    ```
    mkdir _build
    cd _build
    cmake .. -DCMAKE_INSTAL_PREFIX=$HOME/hstcal
    ```

2. Build

    ```
    make
    ```

3. Install

    ```
    make install
    ```

4. Add HSTCAL to `$PATH`

    ```bash
    export PATH=$HOME/hstcal/bin:$PATH
    ```


## Build Targets

To install individual parts of HSTCAL...

```
make install acscte
```

To list available build targets:

```
make help
```

Some common targets include:

Target | Description
-------|------------
acs    | calacs
wf3    | calwf3 and other WFC3-related tools
stis   | calstis


## Debugging

To enable support for debugging symbols use one of the following defines, `cmake .. -DCMAKE_BUILD_TYPE=[RelWithDebInfo|Debug]`
