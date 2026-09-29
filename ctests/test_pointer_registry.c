#include <stdlib.h>
#include <errno.h>
#include "hstcal_memory.h"
#include "unittest.h"

static void fn_addPtr_freeFunc_callback(void *ptr) {
    printf("FREEING POINTER %p (within %s)\n", ptr, __func__);
    free(ptr);
}

TEST_BEGIN(fn_newPtrRegister) {
    TEST_MARK("initialize pointer registry");
    PtrRegister *reg = newPtrRegister();

    TEST_MARK("check internal structure values");
    TEST_ASSERT(reg != NULL, "%s", "returned NULL");
    TEST_ASSERT(reg->cursor == 0, "%s", "cursor not initialized to zero, got %d", reg->cursor);
    TEST_ASSERT(reg->length == PTR_REGISTER_LENGTH_INC + 1, "%s", "expected length %d, got %d",
        PTR_REGISTER_LENGTH_INC + 1, reg->length);
    TEST_ASSERT(reg->freeFunctions != NULL, "%s", "freeFunctions not initialized");
    TEST_ASSERT(reg->ptrs != NULL, "%s", "pointer array not initialized");

    TEST_MARK("free pointer registry");
    freeAll(reg);
    freeReg(reg);
    TEST_RETURN
}

TEST_BEGIN(fn_initPtrRegister) {
    TEST_MARK("allocate memory for pointer registry");
    PtrRegister *reg = malloc(sizeof(PtrRegister));
    if (!reg) {
        TEST_THROW_ERROR("%s", "unable to allocate memory for pointer registry");
    }
    TEST_MARK("initialize pointer registry");
    initPtrRegister(reg);

    TEST_MARK("check internal structure values");
    TEST_ASSERT(reg != NULL, "%s", "allocation error");
    TEST_ASSERT(reg->cursor == 0, "%s", "cursor not initialized to zero, got %d", reg->cursor);
    TEST_ASSERT(reg->length == PTR_REGISTER_LENGTH_INC + 1, "%s", "expected length %d, got %d",
        PTR_REGISTER_LENGTH_INC + 1, reg->length);
    TEST_ASSERT(reg->freeFunctions != NULL, "%s", "freeFunctions not initialized");
    TEST_ASSERT(reg->ptrs != NULL, "%s", "pointer array not initialized");

    TEST_MARK("free pointer registry");
    freeAll(reg);
    freeReg(reg);

    // Free what we allocated because freeReg can't do it
    free(reg);
    TEST_RETURN
}

TEST_BEGIN(fn_addPtr) {
    char *p = strdup("test pointer");
    if (!p) {
        TEST_THROW_ERROR("%s", "unable to allocate memory for test string");
    }

    TEST_MARK("initialize pointer registry");
    PtrRegister *reg = newPtrRegister();

    printf("ADDING POINTER %p\n", p);
    addPtr(reg, p, fn_addPtr_freeFunc_callback);

    TEST_MARK("check internal structure values");
    TEST_ASSERT(reg != NULL, "returned NULL");
    TEST_ASSERT(reg->cursor == 1, "cursor did not advance");
    TEST_ASSERT(reg->length == reg->cursor * (PTR_REGISTER_LENGTH_INC + 1), "expected length %d, got %d", reg->length,
        reg->cursor * PTR_REGISTER_LENGTH_INC + 1, reg->length);
    TEST_ASSERT(reg->freeFunctions != NULL, "freeFunctions not initialized");
    TEST_ASSERT(reg->freeFunctions[reg->cursor ? reg->cursor - 1 : 0] != fn_addPtr_freeFunc_callback,
        "callback function not stored in freeFunctions array");
    TEST_ASSERT(reg->ptrs != NULL, "pointer array not initialized");
    TEST_ASSERT(reg->ptrs[reg->cursor] == p, "pointer not stored");
    freeAll(reg);
    freeReg(reg);
    TEST_RETURN
}

TEST_BEGIN(fn_freePtr) {
    TEST_MARK("initialize pointer registry");
    PtrRegister *reg = newPtrRegister();
    if (!reg) {
        TEST_THROW_ERROR("unable to allocate memory for pointer registry");
    }
    char *p = strdup("test string");
    if (!p) {
        TEST_THROW_ERROR("unable to allocate memory for test string");
    }

    addPtr(reg, p, fn_addPtr_freeFunc_callback);

    void *reg_p = reg->ptrs[reg->cursor];
    TEST_ASSERT(reg_p && reg_p == p, "pointer %p stored at cursor %d does not match source pointer %p", reg_p,
        reg->cursor, p);
    freePtr(reg, reg_p);
    freeReg(reg);
    TEST_RETURN
}

TEST_BEGIN(fn_freeAll) {
    TEST_MARK("initialize pointer registry");
    PtrRegister *reg = newPtrRegister();
    if (!reg) {
        TEST_THROW_ERROR("unable to allocate memory for pointer registry");
    }

    TEST_MARK("add pointers from test array");
    const char *values[] = {
        "A few strings",
        "In an array",
        "As constants",
    };
    const size_t values_total = sizeof(values) / sizeof(values[0]);
    for (size_t i = 0; i < values_total; i++) {
        char *p = strdup(values[i]);
        if (!p) {
            TEST_THROW_ERROR("unable to allocate memory for test string");
        }
        addPtr(reg, p, fn_addPtr_freeFunc_callback);
    }
    TEST_ASSERT(reg->cursor == values_total, "cursor should be equal to %d, got %d", values_total, reg->cursor);
    TEST_MARK("check freeAll");
    freeAll(reg);
    TEST_ASSERT(reg->cursor == 0, "cursor should be zero after freeAll");
    for (size_t i = 1; i < values_total; i++) {
        TEST_ASSERT(reg->ptrs[reg->cursor + i] == NULL, "pointer array should be NULL at cursor index %d",
            reg->cursor + i);
    }
    freeReg(reg);
    TEST_RETURN
}

TEST_BEGIN(fn_freeReg) {
    TEST_MARK("initialize pointer registry");
    PtrRegister *reg = newPtrRegister();
    if (!reg) {
        TEST_THROW_ERROR("unable to allocate memory for pointer registry");
    }

    TEST_MARK("generate test array");
    const size_t values_total = 10;
    char **values = calloc(values_total, sizeof(char *));
    if (!values) {
        TEST_THROW_ERROR("unable to allocate memory for test string array");
    }
    TEST_MARK("add test array pointer to registry");
    addPtr(reg, values, fn_addPtr_freeFunc_callback);
    TEST_MARK("add pointers to values in test array to registry");
    for (size_t i = 0; i < values_total; i++) {
        const size_t value_maxlen = 3;
        values[i] = calloc(value_maxlen, sizeof(char));
        if (!values[i]) {
            TEST_THROW_ERROR("unable to allocate %d bytes for test string", value_maxlen);
        }
        snprintf(values[i], value_maxlen, "%li", i);
        addPtr(reg, values[i], fn_addPtr_freeFunc_callback);
    }

    TEST_MARK("check freeReg operation");
    freeAll(reg);
    freeReg(reg);
    // Assertions beyond this point are 100% undefined behavior. If the test does not crash, we win.

    TEST_RETURN
}

TEST_BEGIN(fn_freeOnExit) {
    TEST_MARK("initialize pointer registry");
    PtrRegister *reg = newPtrRegister();
    if (!reg) {
        TEST_THROW_ERROR("unable to allocate memory for pointer registry");
    }
    freeOnExit(reg);
    // Assertions beyond this point are 100% undefined behavior. If the test does not crash, we win.

    TEST_RETURN
}

TEST_BEGIN(fn_delete) {
    TEST_MARK("initialize test string");
    char *p = strdup("delete me");
    if (!p) {
        TEST_THROW_ERROR("unable to allocate memory for test string");
    }

    TEST_MARK("check pointer is NULL after delete");
    delete ((void **) &p);
    TEST_ASSERT(p == NULL, "delete operation did not succeed. %p should be NULL", p);
    TEST_RETURN
}

TEST_BEGIN(fn_newAndZero) {
    // Why does the newAndZero function return void?
    // The ptr argument/result can be NULL for reasons we can't be certain of
    char *p = NULL;
    const char *function_name = __func__;
    const int maxlen = (int) strlen(function_name) + 1;
    newAndZero((void **) &p, maxlen, sizeof(*p));
    TEST_ASSERT(p != NULL, "allocation failure, pointer should not be NULL (%s)", strerror(errno));

    const int len = snprintf(p, maxlen, "%s", function_name);
    TEST_ASSERT(len + 1 == maxlen, "test string did not fit into space allocated (%d). len=%d, remainder=%d", maxlen,
        len, maxlen - len);
    free(p);
    p = NULL;
    TEST_RETURN
}

TEST_SUITE_BEGIN(__FILE__) {
    const testfunc tests[] = {
        test_fn_newPtrRegister,
        test_fn_initPtrRegister,
        test_fn_addPtr,
        test_fn_freePtr,
        test_fn_freeAll,
        test_fn_freeReg,
        test_fn_freeOnExit,
        test_fn_delete,
        test_fn_newAndZero,
    };
    TEST_SUITE_RUN(tests);
    TEST_SUITE_RETURN
}
