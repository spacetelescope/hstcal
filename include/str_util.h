#ifndef STR_UTIL_INCL
#define STR_UTIL_INCL

#include <stdbool.h>

#define DELIM_NOT_FOUND 0
#define DELIM_FOUND 1
#define DELIM_REJECTED 2

int delim_check(const int ch, const char *accept, const char *reject);
int repchar_s(const char ch, size_t max_ch, char * restrict dest, size_t dest_size);
void upperCase(char * str);
bool isStrInLanguage(const char * str, const char * alphabet);

#endif
