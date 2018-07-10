project(cvos)
add_library(${PROJECT_NAME}
	c_iraf_priv.c
       	irafinit.c
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${HSTCAL_include}
)
