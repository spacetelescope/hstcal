project(wf3_top)
add_subdirectory(lib)

project(wf32d)
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


project(wf3ccd)
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


project(wf3cte)
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


project(wf3ir)
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


project(wf3rej)
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


project(wf3sum)
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


project(calwf3)
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
