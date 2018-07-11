project(acsrej)
add_executable(${PROJECT_NAME}
	acsrej.c
	acsrej_check.c
	acsrej_do.c
	acsrej_init.c
	acsrej_loop.c
	acsrej_sky.c
	cr_history.c
	cr_mode.c
	cr_scaling.c
	mainrej.c
	piksrt.c
	readpar.c
	rej_command.c
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${HSTCAL_include}
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${PROJECT_NAME}/../include
)
target_link_libraries(${PROJECT_NAME}
	PUBLIC acs
)
install(TARGETS ${PROJECT_NAME}
	DESTINATION bin
)

