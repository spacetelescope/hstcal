project(hstio)
add_library(${PROJECT_NAME} SHARED
	hstio.c
	keyword.c
	numeric.c
)
target_link_libraries(${PROJECT_NAME}
	PUBLIC cvos
	PUBLIC ${cfitsio_LDFLAGS}
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${cfitsio_INCLUDE_DIRS})
target_include_directories(${PROJECT_NAME}
	PUBLIC ${HSTCAL_include}
)
target_compile_options(${PROJECT_NAME} PUBLIC ${cfitsio_CFLAGS_OTHER})

