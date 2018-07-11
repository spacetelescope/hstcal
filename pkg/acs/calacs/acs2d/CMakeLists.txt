project(acs2d)
add_executable(${PROJECT_NAME}
	acs2d.c
	do2d.c
	dodark.c
	doflash.c
	doflat.c
	dononlin.c
	dophot.c
	doshad.c
	dospot.c
	get2dflags.c
	get2dsw.c
	getlintab.c
	getspottab.c
	main2d.c
	photmode.c
	sanity2d.c
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${HSTCAL_include}
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${PROJECT_NAME}/../include
)
target_link_libraries(${PROJECT_NAME}
	PUBLIC acs
	PUBLIC hstcalib
)
install(TARGETS ${PROJECT_NAME}
	DESTINATION bin
)

