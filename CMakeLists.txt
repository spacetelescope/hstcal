cmake_minimum_required(VERSION 3.11)
project(HSTCAL)
enable_language(C)

find_package(PkgConfig)
pkg_check_modules(cfitsio cfitsio REQUIRED)

#set(CMAKE_MACOSX_RPATH 1)

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# Keep things looking the way they used to
SET(CMAKE_EXECUTABLE_SUFFIX_C ".e")

include_directories(include)
#add_subdirectory(applib)
add_subdirectory(ctegen2)
add_subdirectory(cvos)
add_subdirectory(hstio)
add_subdirectory(lib)
add_subdirectory(tables)
add_subdirectory(pkg)

