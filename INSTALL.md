# Installation

## Supported Platforms

- Linux >=RHEL6
- OS X >=10.5


## Requirements

- cmake >=3.11 (https://cmake.org)
- cfitsio >=3.430 (https://heasarc.gsfc.nasa.gov/fitsio/fitsio.html)
- gcc >=4.4.7 (https://gcc.gnu.org)
- openmp (http://www.openmp.org/)
- pkg-config (https://www.freedesktop.org/wiki/Software/pkg-config)
- python >=2.7 (https://www.python.org)

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

- If `pkg-config` is not installed, define the path to CFITSIO manually via the `--with-cfitsio` argument:

    ```
    mkdir _build
    cd _build
    cmake .. -DWITH_CFITSIO=$HOME/programs/cfitsio
    ```


### MacOS / OS X

The LLVM/Clang suite provided by Apple XCode is not sufficient to compile HSTCAL. Please install GCC (C and Fortran compilers) either from source, or using a package management system such as Homebrew, MacPorts, Fink, or Conda.


## Compiling HSTCAL

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

3. Install [OPTIONAL]

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
