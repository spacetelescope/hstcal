# UNITTEST

A lightweight single-header unit testing framework for C.

## Features

- Single-header library
- Test suites and test cases
- Optional setup and teardown fixtures
- Assertions
- Skip, fail, and error handling
- Capture and examine test output

## Quickstart

```c
#include "unittest.h"

TEST_UNIT(addition) {
    TEST_ASSERT(1 + 1 == 2, "expected 2");
    TEST_RETURN
}

TEST_UNIT(subtraction) {
    TEST_ASSERT(5 - 2 == 3, "expected 3");
    TEST_RETURN
}

TEST_SUITE("Math Operations") {
    unit_test tests[] = {
        TEST_UNIT_REPR(addition),
        TEST_UNIT_REPR(subtraction),
    };
    TEST_SUITE_RUN(tests)
    TEST_SUITE_RETURN
}
```

### Build

```c
gcc -Wall -Wextra -o test_math_operations tests_math_operations.c
```

### Output

```c
./test_math_operations
UNITTEST v1.0.0
Copyright (C) 2026 Association of Universities for Research in Astronomy (AURA)

 SUITE Math Operations...
  UNIT test__addition... PASSED
  UNIT test__subtraction... PASSED

REPORT Math Operations

Tests ... 2
Pass  ... 2
Fail  ... 0
Error ... 0
Skip  ... 0
```

# Assertions

## Basic Assertions

```c
TEST_ASSERT(x != y, "x should never be equal to y")
```

## Formatted Assertions

```c
int expected = 5;
int result = 6;
TEST_ASSERT(result == expected, "expected %d, got %d", expected, result);
```

### Output

```c
ASSERTION FAILED count == expected BECAUSE expected 5, got 6
```

# Test Status Handling

`TEST_ASSERT` increments an internal error counter whenever the tested condition returns `false`.
`TEST_RETURN` returns `TEST_T_PASS` when the internal counter is zero, or `TEST_T_FAIL` otherwise.

## Pass

```c
TEST_UNIT(zero) {
    TEST_ASSERT(0 == 0, "expected zero!");
    TEST_RETURN  // pass
}
```

## Fail

```c
TEST_UNIT(not_zero) {
    TEST_ASSERT(0 == 1, "expected zero!");
    TEST_RETURN  // fail
}
```

## Error

Use `TEST_THROW_ERROR` when a test cannot continue due to an unrecoverable error.

```c
const char *filename = "datafile";
FILE *fp = fopen(datafile, "rb");
if (!fp) {
    TEST_THROW_ERROR("unable to open %s for reading: %s", filename, strerror(errno));
}
```

## Skip

```c
TEST_THROW_SKIP("Not implemented");
```

To conditionally skip a test, use an `if` statement

```c
if (!use_large_dataset) {
    TEST_THROW_SKIP("Test run not configured for large datasets");
}
```

# Fixtures

Setup and teardown fixtures are defined per-suite and run before, and after each unit test.

## Setup/Teardown

```c
static int setup(void) {
    puts("Inside setup");
    return 0;
}

static int teardown(void) {
    puts("Inside teardown");
    return 0;
}
```

## Register Fixtures

```c
TEST_SUITE("Fixtures") {
    TEST_SUITE_SET_FIXTURE_SETUP(setup);
    TEST_SUITE_SET_FIXTURE_TEARDOWN(teardown);
    unit_test tests[] = {
        TEST_UNIT_REPR(mytest)
    }
    TEST_SUITE_RUN(tests);
    TEST_SUITE_RETURN
}
```

### Output

```
 SUITE Fixtures...
  UNIT test__mytest... PASSED
  Inside setup
  Inside teardown
```

# Marking

Use a `TEST_MARK` to document test execution. This aids debugging unexpected crashes (`SIGSEGV`, etc.).

```c
TEST_UNIT(mytest) {
    TEST_MARK("check strstr");
    TEST_ASSERT(strstr("this substring is present", "substring"), "substring should match");
    TEST_MARK("check strncmp");
    TEST_ASSERT(strncmp(badptr + 3, "string", 6) == 0, "always segfaults")
    TEST_MARK("check strpbrk");
    TEST_ASSERT(strpbrk("this has spaces", " ") != NULL, "string should contain spaces");
    TEST_RETURN
}
```

## Output

```
SUITE mysuite...
  UNIT test__mytest... 
  MARK check strstr
  MARK check strncmp
segmentation fault
```

Breadcrumbs from the crash may be available in the unit test's output log.

```
less -R unitest_log/mysuite/test_mytest.c/test__mytest.log
```

# Testing Captured Output

Use `TEST_REDIRECT_SAVE` to save the test log to a file.

```c
static void a_print_only_function(void) {
    printf("Did it work?\n");
    printf("Yes.");
}

TEST_UNIT(a_print_only_function) {
    a_print_only_function();

    const char *filename = "capture.txt";
    if (TEST_REDIRECT_SAVE(&TEST_LOCAL_REDIRECT, filename)) {
        TEST_THROW_ERROR("unable to save captured output to %s", filename);
    }

    size_t line_count = 0;
    char **lines = TEST_FILE_AS_ARRAY(filename, &line_count);
    if (!lines) {
        TEST_THROW_ERROR("unable to read %s", filename);
    }

    const int expected_lines = 2;
    TEST_ASSERT(line_count == expected_lines, "expected %zu lines, got %zu", expected_lines, line_count);
    if (line_count == expected_lines) {
        // Tests that can only run if we have enough data
        TEST_ASSERT(strcmp(lines[0], "Did it work?\n") == 0, "should have worked");
        TEST_ASSERT(strcmp(lines[1], "Yes.") == 0, "should have worked");
    } else {
        // Otherwise do forensics
        TEST_MARK("debug output");
        size_t last_size = 0;
        for (size_t i = 0; i < line_count; i++) {
            printf("line %zu:\n", i + 1);
            hexdump(lines[i], strlen(lines[i]) + 1);
        }
    }

    TEST_ARRAY_FREE(lines, line_count);
}
```

## Output

### Pass

```
 SUITE test_mytest...
  UNIT test__a_print_only_function... PASSED
```

### Fail

And if a failure is forced by commenting out `printf("Yes.");` in `a_print_only_function`...

```
SUITE test_mytest...
  UNIT test__a_print_only_function... FAILED
       ASSERTION FAILED line_count == expected_lines BECAUSE not enough data in output
       MARK debug output
       line 1:
       00000000  44 69 64 20 69 74 20 77  6f 72 6b 3f 0a 00        Did it work?..
```