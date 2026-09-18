project(wf3_top C)
add_subdirectory(lib)

project(cs0 C)
add_executable(${PROJECT_NAME}
	src/cs0.c
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${HSTCAL_include}
	PUBLIC ${PROJECT_NAME}/../include
)
target_link_libraries(${PROJECT_NAME}
	PUBLIC stis
)
install(TARGETS ${PROJECT_NAME}
	DESTINATION bin
)


project(cs1 C)
add_executable(${PROJECT_NAME}
	src/cs1.c
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${HSTCAL_include}
	PUBLIC ${PROJECT_NAME}/../include
)
target_link_libraries(${PROJECT_NAME}
	PUBLIC stis
)
install(TARGETS ${PROJECT_NAME}
	DESTINATION bin
)


project(cs2 C)
add_executable(${PROJECT_NAME}
	src/cs2.c
	src/cs2_command.c
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${HSTCAL_include}
	PUBLIC ${PROJECT_NAME}/../include
)
target_link_libraries(${PROJECT_NAME}
	PUBLIC stis
)
install(TARGETS ${PROJECT_NAME}
	DESTINATION bin
)


project(cs4 C)
add_executable(${PROJECT_NAME}
	src/cs4.c
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${HSTCAL_include}
	PUBLIC ${PROJECT_NAME}/../include
)
target_link_libraries(${PROJECT_NAME}
	PUBLIC stis
)
install(TARGETS ${PROJECT_NAME}
	DESTINATION bin
)


project(cs6 C)
add_executable(${PROJECT_NAME}
	src/cs6.c
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${HSTCAL_include}
	PUBLIC ${PROJECT_NAME}/../include
)
target_link_libraries(${PROJECT_NAME}
	PUBLIC stis
)
install(TARGETS ${PROJECT_NAME}
	DESTINATION bin
)


project(cs7 C)
add_executable(${PROJECT_NAME}
	src/cs7.c
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${HSTCAL_include}
	PUBLIC ${PROJECT_NAME}/../include
)
target_link_libraries(${PROJECT_NAME}
	PUBLIC stis
)
install(TARGETS ${PROJECT_NAME}
	DESTINATION bin
)


project(cs8 C)
add_executable(${PROJECT_NAME}
	src/cs8.c
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${HSTCAL_include}
	PUBLIC ${PROJECT_NAME}/../include
)
target_link_libraries(${PROJECT_NAME}
	PUBLIC stis
)
install(TARGETS ${PROJECT_NAME}
	DESTINATION bin
)


project(cs11 C)
add_executable(${PROJECT_NAME}
	src/cs11.c
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${HSTCAL_include}
	PUBLIC ${PROJECT_NAME}/../include
)
target_link_libraries(${PROJECT_NAME}
	PUBLIC stis
)
install(TARGETS ${PROJECT_NAME}
	DESTINATION bin
)


project(cs12 C)
add_executable(${PROJECT_NAME}
	src/cs12.c
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${HSTCAL_include}
	PUBLIC ${PROJECT_NAME}/../include
)
target_link_libraries(${PROJECT_NAME}
	PUBLIC stis
)
install(TARGETS ${PROJECT_NAME}
	DESTINATION bin
)
