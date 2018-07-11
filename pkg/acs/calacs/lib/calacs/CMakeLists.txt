project(calacs)
add_executable(${PROJECT_NAME}
	calacs.c
	acsdth.c
	acsinit.c
	acsmain.c
	acsmember.c
	acstable.c
	getflags.c
	getnewref.c
	getsciinfo.c
	refexist.c
	sciflags.c
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${PROJECT_NAME_SOURCE_DIR}
	PUBLIC ${PROJECT_NAME}/../include
	PUBLIC ${HSTCAL_include}
)
target_link_libraries(${PROJECT_NAME}
	PUBLIC acs
	PUBLIC ctegen2
	PUBLIC hstcalib
)
install(TARGETS ${PROJECT_NAME}
	DESTINATION bin
)

