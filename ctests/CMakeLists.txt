add_executable(testphot
    test_getphttab.c
)
add_test(NAME testphot
    COMMAND $<TARGET_FILE:testphot>
)
target_link_libraries(testphot
    PUBLIC hstcalib
)
file(COPY w3m17171j_extrap_imp.fits
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
)

add_executable(test_acscte_rotateamp
    test_acscte_rotateamp.c
)
add_test(NAME test_acscte_rotateamp
    COMMAND $<TARGET_FILE:test_acscte_rotateamp>
)
target_link_libraries(test_acscte_rotateamp
    PUBLIC acs
    PUBLIC hstcalib
)
