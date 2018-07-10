project(hstcalib)
pkg_check_modules(cfitsio cfitsio>=3.2.8 REQUIRED)
add_library(${PROJECT_NAME} STATIC
	hstcal_memory.c
	hstcalversion.c
	trlbuf.c
)
target_link_libraries(${PROJECT_NAME}
	${cfitsio_LIBRARIES}
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${HSTCAL_include}
)
install(TARGETS ${PROJECT_NAME}
	DESTINATION lib
)

