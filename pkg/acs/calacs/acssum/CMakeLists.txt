project(acssum)
add_executable(${PROJECT_NAME}
	acssum.c
	mainsum.c
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

