#ifndef __BOOP_H
#define __BOOP_H
void q(int l, REAL* restrict box, tPart* restrict particles, REAL cutoff, REAL* restrict result);
REAL crystallinity(REAL* restrict box, tPart* restrict particles, REAL cutoff, REAL crystcut);
#endif

