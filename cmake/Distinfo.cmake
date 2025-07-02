function(write_distinfo)
    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/DISTINFO"
            "${GIT_TAG}\n${GIT_BRANCH}\n${GIT_COMMIT}\n"
    )
endfunction()

function(get_distinfo)
    if(EXISTS "${CMAKE_SOURCE_DIR}/DISTINFO")
        file(STRINGS "${CMAKE_SOURCE_DIR}/DISTINFO" DISTINFO)
        set(GIT_TAG DISTINFO[0] PARENT_SCOPE)
        set(GIT_BRANCH DISTINFO[1] PARENT_SCOPE)
        set(GIT_COMMIT DISTINFO[2] PARENT_SCOPE)
    endif()
endfunction()