project(ncarfft)
enable_language(Fortran)
add_library(${PROJECT_NAME}
	ncarfft.f
)
include_directories(${PROJECT_SOURCE_DIR})
install(TARGETS ${PROJECT_NAME}
	DESTINATION lib
)

