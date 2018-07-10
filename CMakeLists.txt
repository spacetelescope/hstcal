cmake_minimum_required(VERSION 3.0)
project(HSTCAL)

find_package(PkgConfig)

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

include_directories(include)
add_subdirectory(applib)
add_subdirectory(cvos)
add_subdirectory(hstio)
add_subdirectory(lib)
add_subdirectory(tables)

#set(hstio_lib_sources
#	hstio/hstio.c
#	hstio/keyword.c
#	hstio/numeric.c
#	)
#
#set(hstio_lib_includes
#	include
#	)
#
#file(GLOB_RECURSE tables_lib_sources
#	tables/*.c
#	)
#
#set(tables_lib_includes
#	include
#	)
#
#
#
#file(GLOB_RECURSE calacs_sources
#	pkg/acs/calacs
#	pkg/acs/calacs/acs*/*.c
#	)
#
#file(GLOB_RECURSE calacs_lib_sources
#	pkg/acs/calacs/lib/*.c
#	)
#
#set(calacs_includes
#	pkg/acs/calacs/include
#	pkg/acs/calacs/acs2d
#	pkg/acs/calacs/acsccd
#	pkg/acs/calacs/acscte
#	pkg/acs/calacs/acsrej
#	pkg/acs/calacs/acssum
#	pkg/acs/calacs/calacs
#	)
#
#set(calacs_lib_includes
#	pkg/acs/calacs/include
#	)
#
#set(calwf3_includes
#	pkg/wfc3/calwf3/include
#	pkg/wfc3/calwf3/wf32d
#	pkg/wfc3/calwf3/wf3ccd
#	pkg/wfc3/calwf3/wf3cte
#	pkg/wfc3/calwf3/wf3ir
#	pkg/wfc3/calwf3/wf3rej
#	pkg/wfc3/calwf3/wf3sum
#	)
#
#file(GLOB_RECURSE calwf3_sources
#	pkg/wfc3/calwf3/calwf3/*.c
#	)
#
#file(GLOB_RECURSE calwf3_lib_sources
#	pkg/wfc3/calwf3/*/*.c
#	)
#
#set(calwf3_lib_includes
#	pkg/wfc3/calwf3/include
#	)
#
#pkg_check_modules(cfitsio cfitsio>=3.2.0 REQUIRED)
#
#include_directories(include)
#
#add_subdirectory(cvos)
##add_library(calacs STATIC ${calacs_lib_sources})
##target_include_directories(calacs PUBLIC ${calacs_lib_includes} ${calacs_includes})
#
#add_library(hstio SHARED ${hstio_lib_sources})
#target_include_directories(hstio PUBLIC ${hstio_lib_includes} ${cfitsio_INCLUDES})
#target_link_libraries(hstio PUBLIC ${cfitsio_LIBRARIES})
#target_link_libraries(hstio PUBLIC tables)
#
#add_library(tables SHARED ${tables_lib_sources})
#target_include_directories(tables PUBLIC ${tables_lib_includes})
#target_link_libraries(tables PUBLIC ${cfitsio_LIBRARIES})
#
#add_library(calwf3 SHARED ${calwf3_lib_sources})
#target_include_directories(calwf3 PUBLIC ${calwf3_lib_includes})
#target_link_libraries(calwf3 PUBLIC ${cfitsio_LIBRARIES})
#target_link_libraries(calwf3 PUBLIC hstio)
#target_link_libraries(calwf3 PUBLIC tables)
#target_link_libraries(calwf3 PUBLIC cvos)
#
#add_executable(calwfc3.e ${calwf3_sources} pkg/wfc3/calwf3/wf3main.c)
#target_link_libraries(calwfc3.e PUBLIC calwf3 cfitsio)
#
#add_executable(wf3ccd.e ${calwf3_sources} pkg/wfc3/calwf3/mainccd.c)
#target_link_libraries(wf3ccd.e PUBLIC calwf3 cfitsio)
#
#add_executable(wf3rej.e ${calwf3_sources} pkg/wfc3/calwf3/mainrej.c)
#target_link_libraries(wf3rej.e PUBLIC calwf3 cfitsio)
#
#add_executable(wf3cte.e ${calwf3_sources} pkg/wfc3/calwf3/maincte.c)
#target_link_libraries(wf3cte.e PUBLIC calwf3 cfitsio)
#
#add_executable(wf3ir.e ${calwf3_sources} pkg/wfc3/calwf3/mainir.c)
#target_link_libraries(wf3ir.e PUBLIC calwf3 cfitsio)
#
#add_executable(wf32d.e ${calwf3_sources} pkg/wfc3/calwf3/main2d.c)
#target_link_libraries(wf32d.e PUBLIC calwf3 cfitsio)
#
#add_executable(wf3sum.e ${calwf3_sources} pkg/wfc3/calwf3/mainsum.c)
#target_link_libraries(wf3sum.e PUBLIC calwf3 cfitsio)
#
#install(TARGETS
#	calwfc3.e
#	wf3ccd.e
#	wf3rej.e
#	wf3cte.e
#	wf3ir.e
#	wf32d.e
#	wf3sum.e
#
#	hstio
#	tables
#
#	calwf3
#
#
#	RUNTIME DESTINATION bin
#	LIBRARY DESTINATION lib)
#
