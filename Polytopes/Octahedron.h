#ifndef OCTAHEDRON_H
#define OCTAHEDRON_H
#include "../include/typedefs.h"
#ifdef _MAIN_
//Definitions for main
#define RC 1.85
static const REAL iscrb_d2=1.100642416298209;
static const REAL cscrb_d2=3.301927248894627;
static const uint TriadV1=1,TriadV2=2; //Which linearly independent Vertices to use. 0,5 for Cube. 1,2 for Octahedron
#else
//Definitions for GJK
static const uint vert_neigh[]={1,2,3,4};

static inline REAL dot_p(const REAL* a, const REAL* b){
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static inline uint poly_max(const REAL* restrict vertices, const REAL* restrict axis){
	//////////////////////////////////////////////////////
	//  Find the maximum projection on the axis using a //
	//	not exactly hill-climbing algorithm. Octahedron	//
	//////////////////////////////////////////////////////
	uint curr=0,next=0;
	REAL p=0.0;
	REAL max=dot_p(vertices,axis);
	
	for(uint i=0;i<4;i++){
		next=vert_neigh[i];
		p=dot_p((vertices+3*next),axis);
		if(p>max){
			max=p;
			curr=next;
		}
	}
	if(curr!=0){
		p=dot_p((vertices+3*5),axis);
		if(p>max){
			curr=5;
		}
	}
	return curr;
}

static inline uint poly_min(const REAL* restrict vertices, const REAL* restrict axis){
	//////////////////////////////////////////////////////
	//  Find the minimum projection on the axis using a //
	//	not exactly hill-climbing algorithm. Octahedron //
	//////////////////////////////////////////////////////
	uint curr=0,next=0;
	REAL p=0.0;
	REAL min=dot_p(vertices,axis);
	
	for(uint i=0;i<4;i++){
		next=vert_neigh[i];
		p=dot_p((vertices+3*next),axis);
		if(p<min){
			min=p;
			curr=next;
		}
	}
	if(curr!=0){
		p=dot_p((vertices+3*5),axis);
		if(p<min){
			curr=5;
		}
	}
	return curr;
}
#endif
#endif
