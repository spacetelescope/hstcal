project(wf3_top C)
add_subdirectory(lib)

project(wf32d C)
add_executable(${PROJECT_NAME}
	src/main2d.c
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${HSTCAL_include}
	PUBLIC ${PROJECT_NAME}/../include
)
target_link_libraries(${PROJECT_NAME}
	PUBLIC wf3
)
install(TARGETS ${PROJECT_NAME}
	DESTINATION bin
)


project(wf3ccd C)
add_executable(${PROJECT_NAME}
	src/mainccd.c
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${HSTCAL_include}
	PUBLIC ${PROJECT_NAME}/../include
)
target_link_libraries(${PROJECT_NAME}
	PUBLIC wf3
)
install(TARGETS ${PROJECT_NAME}
	DESTINATION bin
)


project(wf3cte C)
add_executable(${PROJECT_NAME}
	src/maincte.c
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${HSTCAL_include}
	PUBLIC ${PROJECT_NAME}/../include
)
target_link_libraries(${PROJECT_NAME}
	PUBLIC wf3
	PUBLIC ctegen2
	PUBLIC hstcalib
)
install(TARGETS ${PROJECT_NAME}
	DESTINATION bin
)


project(wf3ir C)
add_executable(${PROJECT_NAME}
	src/mainir.c
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${HSTCAL_include}
	PUBLIC ${PROJECT_NAME}/../include
)
target_link_libraries(${PROJECT_NAME}
	PUBLIC wf3
)
install(TARGETS ${PROJECT_NAME}
	DESTINATION bin
)


project(wf3rej C)
add_executable(${PROJECT_NAME}
	src/mainrej.c
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${HSTCAL_include}
	PUBLIC ${PROJECT_NAME}/../include
)
target_link_libraries(${PROJECT_NAME}
	PUBLIC wf3
)
install(TARGETS ${PROJECT_NAME}
	DESTINATION bin
)


project(wf3sum C)
add_executable(${PROJECT_NAME}
	src/mainsum.c
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${HSTCAL_include}
	PUBLIC ${PROJECT_NAME}/../include
)
target_link_libraries(${PROJECT_NAME}
	PUBLIC wf3
)
install(TARGETS ${PROJECT_NAME}
	DESTINATION bin
)


project(calwf3 C)
add_executable(${PROJECT_NAME}
	src/mainwf3.c
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${HSTCAL_include}
	PUBLIC ${PROJECT_NAME}/../include
)
target_link_libraries(${PROJECT_NAME}
	PUBLIC wf3
)
install(TARGETS ${PROJECT_NAME}
	DESTINATION bin
)
