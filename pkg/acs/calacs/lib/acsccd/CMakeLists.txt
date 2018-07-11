project(acsccd)
add_executable(${PROJECT_NAME}
	acsccd.c
	blev_funcs_postsm4.c
	blevdrift.c
	blevfit.c
	doatod.c
	dobias.c
	doblev.c
	#doblev_v504.c	#... what?
	doccd.c
	dosink.c
	findblev.c
	findover.c
	getacsflag.c
	getccdsw.c
	mainccd.c
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${HSTCAL_include}
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${PROJECT_NAME}/../include
)
target_link_libraries(${PROJECT_NAME}
	PUBLIC acs
#	PUBLIC hstcalib
)
install(TARGETS ${PROJECT_NAME}
	DESTINATION bin
)

