#!/usr/bin/env bash

# Creates a new test unit test suite, and automatically links executable with hstcalib.
#
# $ create_unittest_suite mytest mylib
# Creating test_mytest.c
# Appending suite to CMakeLists.txt
#
# $ tail CMakeLists.txt
# add_executable(test_mytest
#     test_mytest.c
# )
# add_test(NAME test_mytest
#     COMMAND $<TARGET_FILE:test_mytest>
# )
# target_link_libraries(test_mytest
#     PUBLIC hstcalib mylib
# )

# Converts stdin to a string of hex values
to_hex() {
    data_in="$(cat -)"
    echo "${data_in}" | tr -d '\n' | tr -d ' ' | tr -d '\t' | xxd -p -c 1000000
}

# Append the cmake snippet to CMakeLists.txt if not present in the file
append_to_cmakelists() {
    hex_cmake_snippet=$(echo "$cmake_snippet" | to_hex)
    hex_cmakefile=$(cat CMakeLists.txt | to_hex)
    if [[ "$hex_cmakefile" != *"${hex_cmake_snippet}"* ]]; then
        echo "Appending suite to CMakeLists.txt"
        echo "$cmake_snippet" >> CMakeLists.txt
    else
        echo "${suite_name} already in CMakeLists.txt"
    fi
}

# Create a new test. If the test exists, do nothing.
copy_test_template() {
    src="_unittest_template.c"
    dest="${suite_name}.c"
    if [[ -f "$dest" ]]; then
        echo "${dest} already exists" >&2
        return 1
    fi

    echo "Creating ${dest}"
    cp -f "$src" "$dest"
    return 0
}


# main

if (( $# < 1 )); then
  echo "$0 {suite_name} [project_linkage ...]" >&2
  exit 1
fi

suite_name="$1"
shift

project_linkage=""
while (( $# > 0 ))
do
    project_linkage+="$1 "
    shift
done

if ! [[ $suite_name =~ ^test_ ]]; then
    suite_name="test_$suite_name"
fi

cmake_snippet="
add_executable(${suite_name}
    ${suite_name}.c
)
add_test(NAME ${suite_name}
    COMMAND $<TARGET_FILE:${suite_name}>
)
target_link_libraries(${suite_name}
    PUBLIC hstcalib ${project_linkage}
)"

copy_test_template
append_to_cmakelists
