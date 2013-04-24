#ifndef TRIANGULARCUPOLA_SCALED_H
#define TRIANGULARCUPOLA_SCALED_H
#include "../include/typedefs.h"
#ifdef _MAIN_
//Definitions for main
#define RC 1.90344231421
static const REAL iscrb_d2=0.0;
static const REAL cscrb_d2=3.58512379725;
static const uint TriadV1=1,TriadV2=2; //Which linearly independent Vertices to use. 0,5 for Cube. 1,2 for Octahedron
#else
//Definitions for GJK
static const uint vert_neigh[][10]={
 {1, 2, 7, 8, 0, 0, 0, 0, 0, 0},
 {2, 6, 5, 0, 0, 0, 0, 0, 0, 0},
 {1, 4, 3, 0, 0, 0, 0, 0, 0, 0},
 {8, 4, 2, 0, 0, 0, 0, 0, 0, 0},
 {5, 3, 2, 0, 0, 0, 0, 0, 0, 0},
 {4, 6, 1, 0, 0, 0, 0, 0, 0, 0},
 {7, 5, 1, 0, 0, 0, 0, 0, 0, 0},
 {6, 8, 0, 0, 0, 0, 0, 0, 0, 0},
 {3, 7, 0, 0, 0, 0, 0, 0, 0, 0}};

static const uint num_neigh[]={4, 3, 3, 3, 3, 3, 3, 2, 2};

static inline REAL dot_p(REAL *a, REAL *b){
	REAL result=0.0;
	for(uint i=0;i<3;i++){
		result+=a[i]*b[i];
	}
	return result;
}

static inline uint poly_max(REAL *vertices, REAL *axis){
	//////////////////////////////////////////////////////
	//  Find the maximum projection on the axis using a //
	//	hill-climbing algorithm							//
	//////////////////////////////////////////////////////
	uint next=0,last=0,curr=0;
	REAL p=0.0;
	REAL max=dot_p(vertices,axis);
	while(1){
		for(uint i=0;i<num_neigh[curr];i++){
			next=vert_neigh[curr][i];
			if(next!=last){
				p=dot_p((vertices+3*next),axis);
				if(p>max){
					max=p;
					last=curr;
					curr=next;
					break;
				}
			}
			if(i==num_neigh[curr]-1) return curr;
		}
	}
}

static inline uint poly_min(REAL *vertices, REAL *axis){
	//////////////////////////////////////////////////////
	//  Find the maximum projection on the axis using a //
	//	hill-climbing algorithm							//
	//////////////////////////////////////////////////////
	uint next=0,last=0,curr=0;
	REAL p=0.0;
	REAL min=dot_p(vertices,axis);
	while(1){
		for(uint i=0;i<num_neigh[curr];i++){
			next=vert_neigh[curr][i];
			if(next!=last){
				p=dot_p((vertices+3*next),axis);
				if(p<min){
					min=p;
					last=curr;
					curr=next;
					break;
				}
			}
			if(i==num_neigh[curr]-1)return curr;
		}
	}
}
#endif
#endif
