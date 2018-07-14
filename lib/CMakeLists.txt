project(hstcalib C Fortran)
add_library(${PROJECT_NAME} SHARED
	ncarfft.f
	getphttab.c
	hstcal_memory.c
	hstcalversion.c
	trlbuf.c
)
target_link_libraries(${PROJECT_NAME}
	hstio
	tables
	${cfitsio_LDFLAGS}
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${HSTCAL_include}
)
install(TARGETS ${PROJECT_NAME}
	DESTINATION lib
)

