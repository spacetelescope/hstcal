
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


#project(acscteforwardmodel)
#add_executable(${PROJECT_NAME}
#	src/mainacscteforwardmodel.c
#)
#target_include_directories(${PROJECT_NAME}
#	PUBLIC ${HSTCAL_include}
#	PUBLIC ${PROJECT_NAME}/../include
#)
#target_link_libraries(${PROJECT_NAME}
#	PUBLIC acs
#	PUBLIC ctegen2
#	PUBLIC hstcalib
#)
#install(TARGETS ${PROJECT_NAME}
#	DESTINATION bin
#)


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
