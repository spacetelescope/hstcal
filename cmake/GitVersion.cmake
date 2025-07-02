function(get_version_from_git)
    find_package(Git QUIET)
    if (NOT Git_FOUND)
        set(GIT_TAG "unknown" PARENT_SCOPE)
        set(GIT_BRANCH "unknown" PARENT_SCOPE)
        set(GIT_COMMIT "unknown" PARENT_SCOPE)
        message(WARNING "Git not found")
        return()
    endif()

    execute_process(
            COMMAND ${GIT_EXECUTABLE} describe --always --dirty --long --tags
            OUTPUT_VARIABLE GIT_TAG
            OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
            OUTPUT_VARIABLE GIT_BRANCH
            OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
            OUTPUT_VARIABLE GIT_COMMIT
            OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    set(GIT_TAG ${GIT_TAG} PARENT_SCOPE)
    set(GIT_BRANCH ${GIT_BRANCH} PARENT_SCOPE)
    set(GIT_COMMIT ${GIT_COMMIT} PARENT_SCOPE)
endfunction()

