project(cvos)
add_library(${PROJECT_NAME} SHARED
	c_iraf_priv.c
       	irafinit.c
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${HSTCAL_include}
)
set_target_properties(${PROJECT_NAME}
	PROPERTIES
	OUTPUT_NAME hst${PROJECT_NAME}
)
install(TARGETS ${PROJECT_NAME}
	DESTINATION lib
)
