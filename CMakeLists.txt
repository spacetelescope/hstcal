cmake_minimum_required(VERSION 3.11)
project(hstcal)
enable_language(C)
set(CMAKE_C_STANDARD 99)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

option(ENABLE_ASAN "Enable AddressSanitizer" OFF)
option(ENABLE_ASAN_RECOVER "When ASAN reports a problem, don't halt execution" OFF)
option(ENABLE_WARNINGS "Enable compiler warnings" ON)
option(ENABLE_OPENMP "Enable OpenMP" ON)
set(WITH_CFITSIO "" CACHE STRING "Path to cfitsio (if empty pkg-config is used)")
set(WITH_CFITSIO_CFLAGS "" CACHE STRING "CFITSIO compiler flags")
set(WITH_CFITSIO_LDFLAGS "-lcfitsio" CACHE STRING "CFITSIO linker flags")

include(GitVersion)
include(Distinfo)

get_version_from_git()
write_distinfo()
if(NOT EXISTS ${CMAKE_SOURCE_DIR}/.git AND EXISTS ${CMAKE_SOURCE_DIR}/DISTINFO)
	read_distinfo()
endif()

set(PROJECT_VERSION "${GIT_TAG}")
message("HSTCAL Version: ${PROJECT_VERSION}")

# Hack to produce autotools-like distribution archives
add_custom_target(dist
	COMMAND ${CMAKE_COMMAND}
		-D BDIR="${CMAKE_CURRENT_BINARY_DIR}"
		-D PROJECT_NAME="${PROJECT_NAME}"
		-D PROJECT_VERSION="${PROJECT_VERSION}"
		-P ${CMAKE_SOURCE_DIR}/cmake/DistinfoArchive.cmake
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

include(CheckSymbolExists)
check_symbol_exists(snprintf "stdio.h" HAVE_SNPRINTF)
check_symbol_exists(strdup "string.h" HAVE_STRDUP)
check_symbol_exists(INT_MAX "limits.h" HAVE_INT_MAX)
check_symbol_exists(PATH_MAX "limits.h" HAVE_LIMITS__PATH_MAX)
check_symbol_exists(PATH_MAX "sys/limits.h" HAVE_SYS_LIMITS__PATH_MAX)
check_symbol_exists(PATH_MAX "sys/syslimits.h" HAVE_SYS_SYSLIMITS__PATH_MAX)

configure_file(config.h.in config.h @ONLY)

if(NOT HAVE_SNPRINTF OR NOT HAVE_STRDUP OR NOT HAVE_INT_MAX)
	message(FATAL_ERROR "Required symbol is missing!")
endif()

if(ENABLE_ASAN)
	# Report memory leaks, buffer overflows, etc.
	# See: https://clang.llvm.org/docs/AddressSanitizer.html
	add_compile_options(-fno-omit-frame-pointer)
	add_compile_options(-fsanitize=address)
	add_link_options(-fsanitize=address)
	if(ENABLE_ASAN_RECOVER)
		# To use this option set the following in your local environment:
		# ASAN_OPTIONS="halt_on_error=0"
		#
		# ASAN will continue past the first detected issue
		add_compile_options(-fsanitize-recover=address)
		add_link_options(-fsanitize-recover=address)
	endif()
endif()

if(ENABLE_WARNINGS)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")
	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wextra")
endif()

# Initialize pkg-config whether it is used or not
find_package(PkgConfig)

# Figure out where cfitsio should come from
if("${WITH_CFITSIO}" STREQUAL "")
	pkg_check_modules(cfitsio cfitsio REQUIRED)
	set(CMAKE_INSTALL_RPATH "${cfitsio_LIBRARY_DIRS}:${CMAKE_INSTALL_RPATH}")
else()
	set(cfitsio_INCLUDE_DIRS "${WITH_CFITSIO}/include")
	set(cfitsio_LIBRARY_DIRS "${WITH_CFITSIO}/lib")
	set(cfitsio_CFLAGS_OTHER "-I${cfitsio_INCLUDE_DIRS} ${WITH_CFITSIO_CFLAGS}")
	set(cfitsio_LDFLAGS "-L${cfitsio_LIBRARY_DIRS} ${WITH_CFITSIO_LDFLAGS}")
	set(CMAKE_INSTALL_RPATH "${cfitsio_LIBRARY_DIRS}:${CMAKE_INSTALL_RPATH}")
endif()

# Keep things looking the way they used to
set(CMAKE_EXECUTABLE_SUFFIX_C ".e")

# Set up include directories
execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/include)
include_directories("${CMAKE_SOURCE_DIR}/include")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/include")
include_directories("${CMAKE_CURRENT_BINARY_DIR}")

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used by installed binaries
if(APPLE)
	set(CMAKE_MACOSX_RPATH 1)
	set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)
else()
	set(CMAKE_INSTALL_RPATH $ORIGIN/../lib:${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_RPATH})
endif()


set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_PACKAGE_NAME ${PROJECT_NAME})
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
set(CPACK_PACKAGE_VENDOR "Association of Universities for Research in Astronomy (AURA)")
set(CPACK_SYSTEM_NAME ${CMAKE_HOST_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.txt")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
set(CPACK_SOURCE_IGNORE_FILES
	.git.*
	build/
	".*~$"
)
set(CPACK_VERBATIM_VARIABLES YES)
include(CPack)

add_subdirectory(ctegen2)
add_subdirectory(cvos)
add_subdirectory(hstio)
add_subdirectory(lib)
add_subdirectory(tables)
add_subdirectory(pkg)

