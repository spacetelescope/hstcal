project(acscte)
add_executable(${PROJECT_NAME}
	acscte.c
	docte.c
	dopcte-gen2.c
	dopcte.c
	getcteflag.c
	getctesw.c
	maincte.c
	pcte_fixycte.c
	pcte_funcs.c
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${HSTCAL_include}
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${PROJECT_NAME}/../include
)
target_link_libraries(${PROJECT_NAME}
	PUBLIC acs
	PUBLIC ctegen2
	PUBLIC hstcalib
)
install(TARGETS ${PROJECT_NAME}
	DESTINATION bin
)

