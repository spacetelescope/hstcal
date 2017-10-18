#ifndef HSTCAL_MEMORY_INCL
#define HSTCAL_MEMORY_INCL

/* Below is a general use pointer register (book keeping) to allow easy cleanup upon exit (pseudo C++ destructor)
 * Use example:
 *
 * PtrRegister ptrReg;
 * PtrRegister initPtrRegister(&ptrReg); // must always initialize before further use
 *
 * void * array = malloc(size);
 * addPtr(&ptrReg, array, &free); // returns if array == NULL so don't need to pre check - this is also self expanding
 *
 * Bool cancel = False;
 * ...
 * do something
 * cancel = True; //uh oh, something went wrong
 * ...
 * if (cancel)
 * {
 *     freeOnExit(&ptrReg); // This frees itself, i.e. the register array holding the ptrs
 *     return;
 * }
 * ...
 *
 * Alternatively, instead of instantiating an object we can use a pointer along with the functions newPtrRegister(),
 * e.g.
 *
 * PtrRegister * reg = newPtrRegister(); //this adds itself (this pointer) and will free itself by calling
 *                                       // freeReg() or freeOnExit()
 *
 * NOTE: This pattern is considered integral to all use and as such internal failed allocations are asserted
 */

#define PTR_REGISTER_LENGTH_INC 10

typedef void (*FreeFunction)(void*); // Only trivial functions accepted

typedef struct {
    unsigned cursor;
    unsigned length;
    void ** ptrs;
    FreeFunction * freeFunctions;
} PtrRegister;

void * newPtrRegister(); //Allocates a PtrRegister, calls initPtrRegister, registers allocated pointer then returns it
void initPtrRegister(PtrRegister * reg); // initializes members, inc. alloc of registers
void addPtr(PtrRegister * reg, void * ptr, void * freeFunc); // self expanding
void freePtr(PtrRegister * reg, void * ptr); // non contracting
void freeOnExit(PtrRegister * reg); //only calls freeAll() followed by freeOnlyReg()
void freeAll(PtrRegister * reg); //frees all ptrs registered (excluding itself)
void freeReg(PtrRegister * reg); //frees ONLY the registers themselves and NOT the pointers in PtrRegister::ptrs

//Other memory related helper functions
void * newAndZero(void ** ptr, const size_t count, const size_t size);
void delete(void ** ptr);

#endif
