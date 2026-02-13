add_executable(testphot
    test_getphttab.c
)
add_test(NAME testphot
    COMMAND $<TARGET_FILE:testphot>
)
target_link_libraries(testphot
    PUBLIC hstcalib
)
