#ifndef STR_UTIL_INCL
#define STR_UTIL_INCL

#include <stdbool.h>

int repchar_s(const char ch, size_t max_ch, char * restrict dest, size_t dest_size);
void upperCase(char * str);
bool isStrInLanguage(const char * str, const char * alphabet);

#endif
