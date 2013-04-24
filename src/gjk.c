#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "../include/gjk.h"
#include "../include/common.h"
#define CA_S (P*(P-1))/2 //Upper Triangular Matrix Size
#define GJK_CACHE_SIZE CA_S / 50

// Lookup table which tells us at which positions in the simplex array
// to get our points a, b, c, d. I.e. if our bits are 0111 -> 7 -> {0, 1, 2}
static const uchar p_pos[16][3] = 
{
	{0, 0, 0}, {0, 0, 0}, {1, 0, 0}, {0, 1, 0},
	{2, 0, 0}, {0, 2, 0}, {1, 2, 0}, {0, 1, 2},
	{3, 0, 0}, {0, 3, 0}, {1, 3, 0}, {0, 1, 3},
	{2, 3, 0}, {0, 2, 3}, {1, 2, 3}, {0, 0, 0}
};

static const uchar s_pos[] = {0, 0, 1, 0, 2, 0, 0, 0, 3}; //Lookup table for single enabled bit position
//_______________________________^__^_____^___________^


typedef struct{
	int size;
	int capacity;
	REAL* sa;
}GJKCache;

// Here we use a kind of vector to store cached axes
// gjk_keys, tells us in which index of the GJKCache array sa
// to find the cached axis.
static GJKCache gjk_cache;
static int gjk_keys[CA_S];

static inline uint key(uint i, uint j){
	return P*i-((i+1)*(i+2))/2+j;
}


typedef struct{
	REAL p[12];//up to 4 points / 3-Simplex
	uchar bits;
	uchar last_sb;
	uchar size;
}Simplex;

static inline uchar addPoint(Simplex* s){
	uchar b = ~s->bits; //Flip bits
	b &= -b; //Last set (available) bit
	uchar pos = s_pos[b]; //Get the bit position from the lookup table
	s->last_sb = pos;
	s->bits |= b; //Insert the new bit
	s->size++;
	return pos;
}

static inline void removePoint(Simplex *s, int p){
	s->bits ^= (1 << p); //Erase the bit at position p
	s->size--;
}

void init_gjk_cache(void){
	//////////////////////////////////////////
	//  Initialize the cached axes array	//
	//////////////////////////////////////////
	gjk_cache.size = 0;
	gjk_cache.capacity = GJK_CACHE_SIZE;
	gjk_cache.sa = malloc(3*gjk_cache.capacity * sizeof(REAL));
	for(uint i=0;i<CA_S;i++){
		gjk_keys[i] = -1;
	}
}

void free_gjk_cache(void){
	//////////////////////////////////////////
	//  Initialize the cached axes array	//
	//////////////////////////////////////////
	free(gjk_cache.sa);
}

static inline void cross_p(
	const REAL* restrict a,
	const REAL* restrict b,
	REAL* restrict c)
{
	//////////////////////
	//  Cross Product	//
	//////////////////////
	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];
}

static inline void triple_p(
	const REAL* a,
	const REAL* b,
	const REAL* c,
	REAL* result)
{
	//////////////////////////////
	// Triple product (axb)xc 	//
	//////////////////////////////
	REAL cross_temp[3];
	cross_p(a, b, cross_temp);
	cross_p(cross_temp, c, result);
}

static inline void support(
	Simplex* restrict s,
	const REAL* restrict vert_a,
	const REAL* restrict vert_b,
	const REAL* restrict real_dist,
	const REAL* restrict dir)
{
	//////////////////////////////////////////////////////
	//  The support function. This provides us with 	//
	//	the next point for our simplex.					//
	//////////////////////////////////////////////////////
	uint pa,pb;
	pa = poly_max(vert_a, dir);
	pb = poly_min(vert_b, dir);
	
	int pos = addPoint(s);
	
	s->p[3*pos+0] = real_dist[0] + vert_a[3*pa+0] - vert_b[3*pb+0];
	s->p[3*pos+1] = real_dist[1] + vert_a[3*pa+1] - vert_b[3*pb+1];
	s->p[3*pos+2] = real_dist[2] + vert_a[3*pa+2] - vert_b[3*pb+2];
}
	
static bool containsOrigin(
	Simplex* restrict s,
	REAL* restrict dir)
{
	///////////////////////////////////////////////
	//  Check if the origin is contained in the  //
	//	Minkowski sum. 							 //
	///////////////////////////////////////////////
	switch( s->size ){
	case 4:
	{
		/////////////////////////////////////////////////////
		///////////////////* Tetrahedron *///////////////////
		/////////////////////////////////////////////////////
		
		const uchar* pos = p_pos[(s->bits ^ (1 << s->last_sb))];
		
		REAL *a = (s->p + 3*s->last_sb);
		REAL *b = (s->p + 3*pos[0]);
		REAL *c = (s->p + 3*pos[1]);
		REAL *d = (s->p + 3*pos[2]);
		REAL ab[3], ac[3], ad[3];
		REAL abcPerp[3], acdPerp[3], abdPerp[3];
		REAL abPerp1[3], abPerp2[3];
		REAL acPerp1[3], acPerp2[3];
		REAL adPerp1[3], adPerp2[3];
		REAL abxac[3], abxad[3], acxad[3];
		
		for(uint i = 0; i < 3; i++){
			ab[i] = b[i] - a[i];
			ac[i] = c[i] - a[i];
			ad[i] = d[i] - a[i];
		}
		
		////////////////////* Face Cases *///////////////////
		
		
		/* Find the triangle face normals with the correct sign (pointing outward) */
		cross_p(ab, ac, abxac);
		cross_p(ab, ad, abxad);
		cross_p(ac, ad, acxad);
		
		int sign_abc = (dot_p(abxac, ad) > 0.0)? -1: 1;
		int sign_acd = (dot_p(acxad, ab) > 0.0)? -1: 1;
		int sign_abd = (dot_p(abxad, ac) > 0.0)? -1: 1;
		for(uint i = 0; i < 3; i++){
			abcPerp[i] = sign_abc * abxac[i];
			acdPerp[i] = sign_acd * acxad[i];
			abdPerp[i] = sign_abd * abxad[i];
		}
		
		cross_p(acxad, ac, acPerp1);
		cross_p(ad, acxad, adPerp2);
		bool acPerp1Pos = (dot_p(acPerp1, a) > 0.0);
		bool adPerp2Pos = (dot_p(adPerp2, a) > 0.0);
		/* On acd side */
		// The origin should be on acd's side and between the half-spaces defined by ac and ad (normal to acd)
		if((dot_p(acdPerp, a) < 0.0) && !acPerp1Pos && !adPerp2Pos){
			/* Remove point b */
			removePoint(s, pos[0]);
			memcpy(dir, acdPerp, 3*sizeof(*dir));
			break;
		}
		
		cross_p(abxad, ab, abPerp2);
		cross_p(ad, abxad, adPerp1);
		bool abPerp2Pos = (dot_p(abPerp2, a) > 0.0);
		bool adPerp1Pos = (dot_p(adPerp1, a) > 0.0);
		/* On abd side */
		// The origin should be on abd's side and between the half-spaces defined by ab and ad (normal to abd)
		if((dot_p(abdPerp, a) < 0.0) && !abPerp2Pos && !adPerp1Pos){
			/* Remove point c */
			removePoint(s, pos[1]);
			memcpy(dir, abdPerp, 3*sizeof(*dir));
			break;
		}
		
		cross_p(abxac, ab, abPerp1);
		cross_p(ac, abxac, acPerp2);
		bool abPerp1Pos = (dot_p(abPerp1, a) > 0.0);
		bool acPerp2Pos = (dot_p(acPerp2, a) > 0.0);
		/* On abc side */
		// The origin should be on abc's side and between the half-spaces defined by ac and ab (normal to abc)
		if((dot_p(abcPerp,a) < 0.0) && !abPerp1Pos && !acPerp2Pos){
			/* Remove point d */
			removePoint(s, pos[2]);
			memcpy(dir, abcPerp, 3*sizeof(*dir));
			break;
		}
		
		////////////////////* Edge Cases *///////////////////
		
		/* ab Edge case */
		// The origin must be inside the space defined by the intersection
		// of two half-space normal to the adjacent faces abc, abd 
		if(abPerp1Pos && abPerp2Pos){
			triple_p(a, ab, ab, dir);
			removePoint(s, pos[1]);
			removePoint(s, pos[2]);
			break;
		}
		
		
		/* ac Edge case */
		// The origin must be inside the space defined by the intersection
		// of two half-space normal to the adjacent faces abc, acd 
		if(acPerp1Pos && acPerp2Pos){
			triple_p(a, ac, ac, dir);
			removePoint(s, pos[0]);
			removePoint(s, pos[2]);
			break;
		}
		
		/* ad Edge case */
		// The origin must be inside the space defined by the intersection
		// of two half-space normal to the adjacent faces acd, abd 
		if(adPerp1Pos && adPerp2Pos){
			triple_p(a, ad, ad, dir);
			removePoint(s, pos[0]);
			removePoint(s, pos[1]);
			break;
		}
		/* 'else' should only be when the origin is inside the tetrahedron */
		return true;
	}
	case 3:
	{
		/////////////////////////////////////////////////////
		/////////////////////* Triangle *////////////////////
		/////////////////////////////////////////////////////
		
		const uchar* pos = p_pos[(s->bits ^ (1 << s->last_sb))];
		
		REAL *a = (s->p + 3*s->last_sb);
		REAL *b = (s->p + 3*pos[0]);
		REAL *c = (s->p + 3*pos[1]);
		REAL ab[3], ac[3];
		REAL abPerp[3], acPerp[3];
		
		for(uint i = 0; i < 3; i++){
			ab[i] = b[i] - a[i];
			ac[i] = c[i] - a[i];
		}
		
		////////////////////* Edge Cases *///////////////////
		
		REAL abxac[3];
		cross_p(ab, ac, abxac);
		
		/* Origin on the outside of triangle and close to ab */
		cross_p(ab, abxac, abPerp);
		if(dot_p(abPerp, a) < 0.0){
			triple_p(a, ab, ab, dir);
			/* Remove Point c */
			removePoint(s, pos[1]);
			break;
		}
		
		/* Origin on the outside of triangle and close to ac */
		cross_p(abxac, ac, acPerp);
		if(dot_p(acPerp, a) < 0.0){
			triple_p(a, ac, ac, dir);
			/* Remove Point b */
			removePoint(s, pos[0]);
			break;
		}
		
		/////////////////////* Face Case *///////////////////
		
		int sign = (dot_p(abxac, a) > 0.0)? -1: 1;
		for(uint i = 0; i < 3; i++) dir[i] = sign * abxac[i];
		break;
	}
	case 2:
	{
		/////////////////////////////////////////////////////
		///////////////////////* Line *//////////////////////
		/////////////////////////////////////////////////////
		
		const uchar* pos = p_pos[(s->bits ^ (1 << s->last_sb))];
		
		REAL *a = (s->p + 3*s->last_sb);
		REAL *b = (s->p + 3*pos[0]);
		REAL ab[3];
		for(uint i = 0; i < 3; i++) ab[i] = b[i] - a[i];
		triple_p(a, ab, ab, dir);
		break;
	}
	case 1:
	{
		/////////////////////////////////////////////////////
		///////////////////////* Point */////////////////////
		/////////////////////////////////////////////////////
		REAL *a = (s->p + 3*s->last_sb);
		for(uint i = 0; i < 3; i++) dir[i] = -a[i]; //Take direction passing through origin
		break;
	}
	default: break;
	}
	return false;
}

bool gjk_overlap(
	const tPart* restrict a,
	const tPart* restrict b,
	const REAL* restrict real_distance)
{
	
	uint indi, indj, hkey;
	
	if(a->index > b->index){ //Upper Triangular Matrix
		indj = a->index;
		indi = b->index;
	}
	else{
		indi = a->index;
		indj = b->index;
	}
	hkey = key(indi,indj);
	/////////////////////////GJK Starts Here//////////////////////////////////
	REAL dir[3] = {1.0, 0.0, 0.0};
	Simplex S;
	S.size = 0;
	S.bits = 0;
	if(gjk_keys[hkey] != -1) memcpy(dir, (gjk_cache.sa + 3*gjk_keys[hkey]), 3*sizeof(*dir));
	
	uint fail_safe = 0;
	
	while(fail_safe < 20){
	
		fail_safe++;
		
		support(&S, a->vertices, b->vertices, real_distance, dir);
		
		if(dot_p((S.p+3*(S.last_sb)), dir) < 0.0){
			if(gjk_keys[hkey] == -1){
				gjk_keys[hkey] = gjk_cache.size;
				gjk_cache.size++;
				if(gjk_cache.size >= gjk_cache.capacity){
					gjk_cache.capacity *= 2;
					gjk_cache.sa = realloc(gjk_cache.sa, 3*gjk_cache.capacity*sizeof(REAL));
				}
			}
			memcpy((gjk_cache.sa + 3*gjk_keys[hkey]), dir, 3*sizeof(*dir));
			return false;
		}
		else if(containsOrigin(&S, dir)){
            return true;
        }
	}
	printf("Encountered error in GJK: Infinite Loop.\n Direction (%f, %f, %f)\n", dir[0], dir[1], dir[2]);
	return false;
}
