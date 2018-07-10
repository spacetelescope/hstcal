project(hstio)
pkg_check_modules(cfitsio cfitsio>=3.2.8 REQUIRED)
add_library(${PROJECT_NAME}
	hstio.c
	keyword.c
	numeric.c
)
target_link_libraries(${PROJECT_NAME} ${cfitsio_LIBRARIES})
target_include_directories(${PROJECT_NAME}
	PUBLIC ${cfitsio_INCLUDE_DIRS})
target_include_directories(${PROJECT_NAME}
	PUBLIC ${HSTCAL_include}
)
target_compile_options(${PROJECT_NAME} PUBLIC ${cfitsio_CFLAGS_OTHER})

