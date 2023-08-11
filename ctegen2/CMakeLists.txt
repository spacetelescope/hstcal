project(ctegen2)
find_package(OpenMP COMPONENTS C)

# Circular dependency on ACS, so this MUST be STATIC
add_library(${PROJECT_NAME} STATIC
	ctegen2.c
	ctehelpers.c
)
target_link_libraries(${PROJECT_NAME}
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

target_compile_options(${PROJECT_NAME}
	PUBLIC "-fPIC"
)

if (OpenMP_FOUND)
	target_link_libraries(${PROJECT_NAME}
		PUBLIC ${OpenMP_C_LIB_NAMES}
	)
	target_compile_options(${PROJECT_NAME}
		PUBLIC ${OpenMP_C_FLAGS}
	)
endif()

set_target_properties(${PROJECT_NAME}
	PROPERTIES
	OUTPUT_NAME hst${PROJECT_NAME}
)
install(TARGETS ${PROJECT_NAME}
	DESTINATION lib
)