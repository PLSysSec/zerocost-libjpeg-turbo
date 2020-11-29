#pragma once

#include <stdlib.h>

typedef int jmp_buf;

inline int setjmp (jmp_buf env) { return 0; }

void longjmp (jmp_buf env, int val) { abort(); }