#include <stdlib.h>
#include <assert.h>

#include "hstio.h"
#include "hstcal_memory.h"

void * newPtrRegister()
{
    void * this = malloc(sizeof(PtrRegister));
    initPtrRegister(this);
    addPtr(this, this, &free); // Note: freeFunctions[0] is ignored anyhow
    return this;
}
void initPtrRegister(PtrRegister * reg)
{
    reg->cursor = 0; //points to last ptr NOT next slot
    reg->length = PTR_REGISTER_LENGTH_INC+1; //+1 to special case reg->ptrs[0] for 'this' pointer only
    reg->ptrs = malloc(reg->length*sizeof(*reg->ptrs));
    assert(reg->ptrs);
    reg->freeFunctions = malloc(reg->length*sizeof(*reg->freeFunctions));
    if (!reg->freeFunctions)
    {
        free(reg->ptrs);
        assert(0);
    }
    reg->ptrs[0] = NULL; //initialize to check against later
}
void addPtr(PtrRegister * reg, void * ptr, void * freeFunc)
{
    if (!reg || !ptr || !freeFunc)
        return;

    //check ptr isn't already registered? - go on then.
    {int i;
    for (i = reg->cursor; i >= 0 ; --i)// i >= 0 prevents adding self again
    {
        if (reg->ptrs[i] == ptr)
            return;
    }}

    if (ptr == reg)
    {
        reg->ptrs[0] = ptr;
        reg->freeFunctions[0] = freeFunc;
        return; //don't inc reg->cursor
    }

    if (++reg->cursor >= reg->length)
    {
        reg->length += PTR_REGISTER_LENGTH_INC;
        void * tmpPtr = realloc(reg->ptrs, reg->length*sizeof(*reg->ptrs));
        if (!tmpPtr)
        {
            freeOnExit(reg); // Note: It is ok that reg->length != length(reg->ptrs)
            assert(0);
        }
        reg->ptrs = tmpPtr;
        tmpPtr = realloc(reg->freeFunctions, reg->length*sizeof(*reg->freeFunctions));
        if (!tmpPtr)
        {
            freeOnExit(reg); // Note: It is ok that reg->length != length(reg->freeFunctions)
            assert(0);
        }
        reg->freeFunctions = tmpPtr;
    }
    reg->ptrs[reg->cursor] = ptr;
    reg->freeFunctions[reg->cursor] = freeFunc;
}
void freePtr(PtrRegister * reg, void * ptr)
{
    //Can't be used to free itself, use freeReg(), use of i > 0 in below for is reason.
    if (!reg || !ptr || !reg->cursor)
        return;

    int i;
    Bool found = False;
    for (i = reg->cursor; i > 0 ; --i)
    {
        if (reg->ptrs[i] == ptr)
        {
            found = True;
            break;
        }
    }
    if (!found)
    {
        freeOnExit(reg); // Note: Whilst the point of this outer IF accounts for direct calls to freePtr() with unregistered ptrs,
                         // it is also called internally from freeOnExit(). A bug in the register code is liable to create an infinite
                         // recursion segfault due this call present here.
        assert(0); //internal error: the ptr trying to be freed was not added to the register
    }

    //call function to free ptr
    reg->freeFunctions[i](ptr);

    if (i == reg->cursor)
    {
        reg->ptrs[i] = NULL;
        reg->freeFunctions[i] = NULL;
    }
    else
    {
        //move last one into gap to close - not a stack so who cares
        reg->ptrs[i] = reg->ptrs[reg->cursor];
        reg->ptrs[reg->cursor] = NULL;
        reg->freeFunctions[i] = reg->freeFunctions[reg->cursor];
        reg->freeFunctions[reg->cursor] = NULL;
    }
    --reg->cursor;
}
void freeAll(PtrRegister * reg)
{
    if (!reg || reg->length == 0 || reg->cursor == 0)
        return;

    while (reg->cursor > 0) //don't free 'this' pointer
        freePtr(reg, reg->ptrs[reg->cursor]);
}
void freeReg(PtrRegister * reg)
{
    /* THIS SHOULD NEVER CALL freeALL()
     * This is designed to be used when allocating multiple persistent memory allocations,
     * registering each allocation in turn. If one allocation fails freeOnExit() can be
     * called to free prior successful allocations and if all allocations are successful
     * this function can be called to free the registers without freeing the actual
     * pointers just allocated.
     */

    if (!reg || reg->length == 0)
        return;

    void * this = reg->ptrs[0];
    // free registers
    free(reg->ptrs);
    reg->ptrs = NULL;
    free(reg->freeFunctions);
    reg->freeFunctions = NULL;

    if (this)
    {
        free(this);
        return;
    }

    reg->cursor = 0;
    reg->length = 0;
}
void freeOnExit(PtrRegister * reg)
{
    //free everything registered
    freeAll(reg);
    //free registers
    freeReg(reg);
}

void delete(void ** ptr)
{
    if (!ptr)
        return;
    if (*ptr)
    {
        free(*ptr);
        *ptr = NULL;
    }
}

void * newAndZero(void ** ptr, size_t count, size_t size)
{
    if (!ptr)
        return NULL;
    delete(ptr);
    *ptr = calloc(count, size);
    return *ptr;
}
