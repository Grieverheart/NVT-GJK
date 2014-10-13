#ifndef QUATERNION_H
#define QUATERNION_H
#include "typedefs.h"

typedef struct{
	REAL el[4];
}quat_t;
void quat_rot(quat_t a, REAL *b, uint nVertices);
#endif
