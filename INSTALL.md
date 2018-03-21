# Installation

HSTCAL uses the WAF build system. A precompiled copy of WAF has been included at the top-level of this project for your convenience.


## Supported Platforms

- Linux >=RHEL6
- OS X >=10.5


## Requirements

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
    ./waf configure --with-cfitsio=$HOME/programs/cfitsio
    ```


### MacOS / OS X

The LLVM/Clang suite provided by Apple XCode is not sufficient to compile HSTCAL. Please install GCC (C and Fortran compilers) either from source, or using a package management system such as Homebrew, MacPorts, Fink, or Conda.


## Compiling HSTCAL

_For a complete listing of useful configuration and build options, run_ `./waf --help`

1. Configure

    ```
    ./waf configure --prefix=$HOME/hstcal
    ```

2. Build

    ```
    ./waf build
    ```

3. Test [OPTIONAL]

    ```
    ./waf test
    ```

4. Install [OPTIONAL]

    ```
    ./waf install
    ```

5. Add HSTCAL to `$PATH`

    ```bash
    export PATH=$HOME/hstcal/bin:$PATH
    ```


## Build Targets

To install individual parts of HSTCAL, use the `--targets` option. For example, to install the `acs` and `lib` targets, do the following:

```
./waf install --targets=acs,lib
```

To list available build targets:

```
./waf configure list
```

Some common targets include:

Target | Description
-------|------------
acs    | calacs
wf3    | calwf3 and other WFC3-related tools
stis   | calstis
lib    | static libraries and header files for the included libraries
test   | self-test executables


## Debugging

To enable support for debugging symbols run, `./waf configure --debug`


## build.cfg file

Shell arguments normally passed to `./waf configure` may be issued via `build.cfg`.  This configuration file may contain any arguments accepted by `./waf configure`.  An example is given in the `build.cfg.example` file.


# Notes for developers


## Using headers


C has two forms of `#include` syntax, `#include "foo.h"` (i.e. `include/foo.h`) for local include files, and `#include <foo.h>` (i.e. `/usr/include/foo.h`) for system include files.  WAF does not track changes to system include files, so the first syntax should be used whenever including files within the project.
