#include "unittest.h"
#include "hstio.h"

// Load symbols that should be exported but aren't
void error(HSTIOError e, char *str);
void ioerr(HSTIOError e, IODescPtr x_, int status);

static void hstio_error_handler() {
#ifdef DEBUG
    TEST_MSG(stderr, NULL, __func__, "called with status %d:\nmessage: '%s'\nmessage length: %zu\n", hstio_err(),
        hstio_errmsg(), strlen(hstio_errmsg()));
#endif
}
static const int handlers_max = 2; // up to 32 handlers, see hstio.h

TEST_BEGIN(macro_Pix) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(macro_PixColumnMajor) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(macro_PPixColumnMajor) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(macro_DQPix) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(macro_DQSetPix) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_copyDataSection) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_hstio_err) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_hstio_errmsg) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_push_hstioerr) {
    TEST_MARK("check state after pushing handler onto stack");
    for (int i = 0; i < handlers_max; i++) {
        const int stack_level = push_hstioerr(hstio_error_handler);
        TEST_ASSERT(stack_level == i + 1, "push failed, got %d", stack_level);
    }
    TEST_RETURN
}

TEST_BEGIN(fn_pop_hstioerr) {
    TEST_MARK("check state after poppping handler from stack");
    for (int i = handlers_max; i >= 0; i--) {
        const int stack_level = pop_hstioerr();
        TEST_ASSERT(stack_level == i - 1, "pop failed, got %d", stack_level);
    }
    TEST_RETURN
}

TEST_BEGIN(fn_clear_hstioerr) {
    TEST_MARK("check error state after clear");
    clear_hstioerr();
    const int err = hstio_err();
    const char *msg = hstio_errmsg();
    TEST_MARK("check error code is zero");
    TEST_ASSERT(err == 0, "expected error code to be zero, got %d\n", err);
    TEST_MARK("check error message is zero-length");
    TEST_ASSERT(strlen(msg) == 0, "expected message length to be zero, got %d\n", strlen(msg));
    TEST_RETURN
}

TEST_BEGIN(fn_error) {
    size_t error_len = 0;
    push_hstioerr(hstio_error_handler);

    TEST_MARK("check contents of error_msg buffer");
    for (size_t i = 0; i <= BADREMOVE; i++) {
        error(i, "canary");
        const char *error_msg = hstio_errmsg();
        TEST_ASSERT(strstr(error_msg, "canary") != NULL, "Message does not contain expected value! error_msg = '%s'\n",
            error_msg);
    }

    TEST_MARK("check unhandled error");
    error(9999, "this is an unhandled error");
    const char *error_msg = hstio_errmsg();
    TEST_ASSERT(strstr(error_msg, "HSTIOError") != NULL, "Failed to trap unhandled error code! error_msg = '%s'\n",
        error_msg);

    TEST_MARK("check very long message truncation");
    char huge_error_message[4096] = {0};
    memset(huge_error_message, '?', sizeof(huge_error_message) - 1);
    error(BADNAME, huge_error_message);
    error_msg = hstio_errmsg();
    error_len = strlen(error_msg);
    const char *expected = "Keyword name";
    TEST_ASSERT(strstr(error_msg, expected) != NULL,
        "Message does not contain expected value, '%s'! error_msg = '%s'\n", expected, error_msg);

    TEST_MARK("HSTOK with an empty error string produces an empty string");
    error(HSTOK, "");
    error_msg = hstio_errmsg();
    error_len = strlen(error_msg);
    TEST_ASSERT(error_len == 0, "Message should be empty! error_msg = '%s'\n", error_msg);

    pop_hstioerr();
    TEST_RETURN
}

TEST_BEGIN(fn_openSingleGroupLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_closeSingleGroupLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_fcloseNull) {
    FILE *fp = fopen("/dev/null", "w");
    if (!fp) {
        //TEST_("%s", "Error opening /dev/null");
        TEST_FORCE_ERROR;
    }
    TEST_MARK("check closing a file handle");
    TEST_ASSERT(fcloseNull(fp) == 0, "%s", "stream did not close\n");

    TEST_MARK("check passing NULL file handle returns zero");
    TEST_ASSERT(fcloseNull(NULL) == 0, "%s", "function returned non-zero on NULL file handle\n");
    TEST_RETURN
}

TEST_BEGIN(fn_fcloseWithStatus) {
    FILE *fp = fopen("/dev/null", "w");
    if (!fp) {
        TEST_THROW_ERROR("%s", "Error opening /dev/null");
    }
    TEST_MARK("check closing a file handle");
    TEST_ASSERT(fcloseWithStatus(&fp) == 0, "%s", "stream did not close\n");

    TEST_MARK("check file handle is NULL");
    TEST_ASSERT(fp == NULL, "%s", "stream is not NULL\n");

    // fcloseWithStatus returns IO_ERROR when fclose() fails. Intentionally
    // forcing fclose() to fail is difficult. Accessing a stream after you've
    // called fclose() on it is purely undefined behavior.
    TEST_RETURN
}

TEST_BEGIN(fn_ckNewFile) {
    TEST_MARK("generating test file");
    const char *filename = "testnewfile";
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        TEST_THROW_ERROR("Error opening %s\n", filename);
    }
    fprintf(fp, "hello world from %s\n", filename);
    fclose(fp);

    int exists = 0;
    exists = ckNewFile((char *) filename);
    TEST_MARK("check non-zero return value if file exists");
    TEST_ASSERT(exists == 1, "existence check failed, returned %d\n", exists);

    TEST_MARK("removing test file");
    remove(filename);
    exists = ckNewFile((char *) filename);
    TEST_MARK("check zero return value if file does not exist");
    TEST_ASSERT(exists == 0, "existence check failed, returned %d\n", exists);
    TEST_RETURN
}

TEST_BEGIN(fn_getSci) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putSci) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getErr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putErr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getDQ) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putDQ) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getSmpl) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putSmpl) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getIntg) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putIntg) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putSciSect) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putErrSect) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putDQSect) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putSmplSect) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putIntgSect) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getSciHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getErrHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getDQHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getSciLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getErrLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getDQLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getFloatHD) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putFloatHD) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getShortHD) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putShortHD) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putFloatHDSect) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putShortHDSect) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getFloatHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getShortHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getSingleGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getSingleGroupLine) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putSingleGroupHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putSingleGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putSingleGroupSect) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getSingleNicmosGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putSingleNicmosGroupHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putSingleNicmosGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putSingleNicmosGroupSect) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getMultiGroupHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getMultiGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putMultiGroupHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putMultiGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getMultiNicmosGroupHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_getMultiNicmosGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putMultiNicmosGroupHdr) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}

TEST_BEGIN(fn_putMultiNicmosGroup) {
    TEST_THROW_SKIP("Stub");
    TEST_RETURN
}


TEST_SUITE_BEGIN(__FILE__) {
    const testfunc tests[] = {
        test_macro_Pix,
        test_macro_PixColumnMajor,
        test_macro_PPixColumnMajor,
        test_macro_DQPix,
        test_macro_DQSetPix,
        test_fn_copyDataSection,
        test_fn_hstio_errmsg,
        test_fn_hstio_err,
        test_fn_push_hstioerr,
        test_fn_pop_hstioerr,
        test_fn_clear_hstioerr,
        test_fn_error,
        test_fn_openSingleGroupLine,
        test_fn_closeSingleGroupLine,
        test_fn_fcloseNull,
        test_fn_fcloseWithStatus,
        test_fn_ckNewFile,
        test_fn_getSci,
        test_fn_putSci,
        test_fn_getErr,
        test_fn_putErr,
        test_fn_getDQ,
        test_fn_putDQ,
        test_fn_getSmpl,
        test_fn_putSmpl,
        test_fn_getIntg,
        test_fn_putIntg,
        test_fn_putSciSect,
        test_fn_putErrSect,
        test_fn_putDQSect,
        test_fn_putSmplSect,
        test_fn_putIntgSect,
        test_fn_getSciHdr,
        test_fn_getErrHdr,
        test_fn_getDQHdr,
        test_fn_getSciLine,
        test_fn_getErrLine,
        test_fn_getDQLine,
        test_fn_getFloatHD,
        test_fn_putFloatHD,
        test_fn_getShortHD,
        test_fn_putShortHD,
        test_fn_putFloatHDSect,
        test_fn_putShortHDSect,
        test_fn_getFloatHdr,
        test_fn_getShortHdr,
        test_fn_getSingleGroup,
        test_fn_getSingleGroupLine,
        test_fn_putSingleGroupHdr,
        test_fn_putSingleGroup,
        test_fn_putSingleGroupSect,
        test_fn_getSingleNicmosGroup,
        test_fn_putSingleNicmosGroupHdr,
        test_fn_putSingleNicmosGroup,
        test_fn_putSingleNicmosGroupSect,
        test_fn_getMultiGroupHdr,
        test_fn_getMultiGroup,
        test_fn_putMultiGroupHdr,
        test_fn_putMultiGroup,
        test_fn_getMultiNicmosGroupHdr,
        test_fn_getMultiNicmosGroup,
        test_fn_putMultiNicmosGroupHdr,
        test_fn_putMultiNicmosGroup,
    };
    TEST_SUITE_RUN(tests);
    TEST_SUITE_RETURN
}
