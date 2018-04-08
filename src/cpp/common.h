#ifndef WASMSOLVER_COMMON_H
#define WASMSOLVER_COMMON_H

#ifdef EMSCRIPTEN
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

#endif