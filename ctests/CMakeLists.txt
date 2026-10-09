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

add_executable(test_hstio
    test_hstio.c
)
add_test(NAME test_hstio
    COMMAND $<TARGET_FILE:test_hstio>
)
target_link_libraries(test_hstio
    PUBLIC hstio
)

add_executable(test_pointer_registry
        test_pointer_registry.c
)
add_test(NAME test_pointer_registry
        COMMAND $<TARGET_FILE:test_pointer_registry>
)
target_link_libraries(test_pointer_registry
        PUBLIC hstcalib
)

add_executable(test_str_util
        test_str_util.c
)
add_test(NAME test_str_util
        COMMAND $<TARGET_FILE:test_str_util>
)
target_link_libraries(test_str_util
        PUBLIC hstcalib
)

add_executable(test_trlbuf
        test_trlbuf.c
)
add_test(NAME test_trlbuf
        COMMAND $<TARGET_FILE:test_trlbuf>
)
target_link_libraries(test_trlbuf
        PUBLIC hstcalib
)


add_executable(test_tables
    test_tables.c
)
add_test(NAME test_tables
    COMMAND $<TARGET_FILE:test_tables>
)
target_link_libraries(test_tables
    PUBLIC hstcalib tables 
)

add_executable(test_timestamp
    test_timestamp.c
)
add_test(NAME test_timestamp
    COMMAND $<TARGET_FILE:test_timestamp>
)
target_link_libraries(test_timestamp
    PUBLIC hstcalib 
)
