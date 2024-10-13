
project(acs_top)
add_subdirectory(lib)

project(acs2d)
add_executable(${PROJECT_NAME}
	src/main2d.c
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${HSTCAL_include}
	PUBLIC ${PROJECT_NAME}/../include
)
target_link_libraries(${PROJECT_NAME}
	PUBLIC acs
)
install(TARGETS ${PROJECT_NAME}
	DESTINATION bin
)


project(acsccd)
add_executable(${PROJECT_NAME}
	src/mainccd.c
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${HSTCAL_include}
	PUBLIC ${PROJECT_NAME}/../include
)
target_link_libraries(${PROJECT_NAME}
	PUBLIC acs
)
install(TARGETS ${PROJECT_NAME}
	DESTINATION bin
)


project(acscte)
add_executable(${PROJECT_NAME}
	src/maincte.c
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${HSTCAL_include}
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

# Handle acscteforwardmodel.e link creation
add_custom_target(${PROJECT_NAME}forwardmodel ALL
	COMMAND ${CMAKE_COMMAND} -E create_hardlink ${PROJECT_NAME}.e ${PROJECT_NAME}forwardmodel.e
)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}forwardmodel.e
	PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
	DESTINATION bin
)


project(acsrej)
add_executable(${PROJECT_NAME}
	src/mainrej.c
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${HSTCAL_include}
	PUBLIC ${PROJECT_NAME}/../include
)
target_link_libraries(${PROJECT_NAME}
	PUBLIC acs
	PUBLIC hstcalib
)
install(TARGETS ${PROJECT_NAME}
	DESTINATION bin
)


project(acssum)
add_executable(${PROJECT_NAME}
	src/mainsum.c
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${HSTCAL_include}
	PUBLIC ${PROJECT_NAME}/../include
)
target_link_libraries(${PROJECT_NAME}
	PUBLIC acs
)
install(TARGETS ${PROJECT_NAME}
	DESTINATION bin
)


project(calacs)
add_executable(${PROJECT_NAME}
	src/mainacs.c
)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${HSTCAL_include}
	PUBLIC ${PROJECT_NAME}/../include
)
target_link_libraries(${PROJECT_NAME}
	PUBLIC acs
)
install(TARGETS ${PROJECT_NAME}
	DESTINATION bin
)
