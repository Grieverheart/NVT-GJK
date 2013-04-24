#ifndef GJK_H
#define GJK_H
#include <stdbool.h>
#include "typedefs.h"
bool gjk_overlap(const tPart* restrict a, const tPart* restrict b, const REAL* restrict real_distance);
void init_gjk_cache(void);
void free_gjk_cache(void);
#endif
