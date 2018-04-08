#include "api.h"
#include "common.h"

EXPORT void *memory_alloc(length_t length) { return new unsigned char[length]; }
EXPORT void memory_free(float *pointer) { delete[](unsigned char *) pointer; }