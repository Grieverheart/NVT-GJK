#ifndef NPT_TYPEDEF_H
#define NPT_TYPEDEF_H
#define REAL double
typedef unsigned int uint;
typedef unsigned short ushort;
typedef unsigned char uchar;
typedef struct{
	REAL center[3];
	REAL real_center[3];
	REAL *vertices;
	ushort index;
}tPart;
#endif
