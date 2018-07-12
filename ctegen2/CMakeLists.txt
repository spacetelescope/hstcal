project(ctegen2)
find_package(OpenMP COMPONENTS C)

add_library(${PROJECT_NAME} STATIC
	ctegen2.c
	ctehelpers.c
)
target_link_libraries(${PROJECT_NAME}
	#PUBLIC acs
	PUBLIC hstio
	PUBLIC tables
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${PROJECT_NAME_SOURCE_DIR}
	PUBLIC ${HSTCAL_include}
	PUBLIC ../pkg/acs/include
	PUBLIC ../pkg/wfc3/include
	PUBLIC ../pkg/stis/include
)

if (OpenMP_FOUND)
	target_link_libraries(${PROJECT_NAME}
		PUBLIC ${OpenMP_C_LIB_NAMES}
	)
	target_compile_options(${PROJECT_NAME}
		PUBLIC ${OpenMP_C_FLAGS}
	)
	target_compile_definitions(${PROJECT_NAME}
		PUBLIC HAVE_OPENMP=1
	)
endif()

install(TARGETS ${PROJECT_NAME}
	DESTINATION lib
)
