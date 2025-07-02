configure_file(${BDIR}/DISTINFO ${CMAKE_SOURCE_DIR}/DISTINFO COPYONLY)
execute_process(
    COMMAND git archive
        -v
        --format=tar.gz
        --output=${BDIR}/${PROJECT_NAME}-${PROJECT_VERSION}-src.tar.gz
        --prefix=${PROJECT_NAME}-${PROJECT_VERSION}/
        --add-file=${BDIR}/DISTINFO
        HEAD .
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
file(REMOVE ${CMAKE_SOURCE_DIR}/DISTINFO)