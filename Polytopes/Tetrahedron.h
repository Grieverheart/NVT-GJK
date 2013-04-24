#ifndef TETRAHEDRON_SCALED_H
#define TETRAHEDRON_SCALED_H
#include "../include/typedefs.h"
#ifdef _MAIN_
//Definitions for main
#define RC 2.6
static const double iscrb_d2=0.693361274349;
static const double cscrb_d2=6.24025146917;
static const uint TriadV1=1,TriadV2=2; //Which linearly independent Vertices to use. 0,5 for Cube. 1,2 for Octahedron
#else
//Definitions for GJK
static const uint vert_neigh[]={1, 2, 3};

static inline double dot_p(double *a, double *b){
	double result=0.0;
	for(uint i=0;i<3;i++){
		result+=a[i]*b[i];
	}
	return result;
}

static inline uint poly_max(double *vertices, double *axis){
	//////////////////////////////////////////////////////
	//  Find the maximum projection on the axis using a //
	//	not exactly hill-climbing algorithm. Octahedron	//
	//////////////////////////////////////////////////////
	uint curr=0,next=0;
	double p=0.0;
	double max=dot_p(vertices,axis);
	
	for(uint i=0;i<3;i++){
		next=vert_neigh[i];
		p=dot_p((vertices+3*next),axis);
		if(p>max){
			max=p;
			curr=next;
		}
	}
	return curr;
}

static inline uint poly_min(double *vertices, double *axis){
	//////////////////////////////////////////////////////
	//  Find the minimum projection on the axis using a //
	//	not exactly hill-climbing algorithm. Octahedron //
	//////////////////////////////////////////////////////
	uint curr=0,next=0;
	double p=0.0;
	double min=dot_p(vertices,axis);
	
	for(uint i=0;i<3;i++){
		next=vert_neigh[i];
		p=dot_p((vertices+3*next),axis);
		if(p<min){
			min=p;
			curr=next;
		}
	}
	return curr;
}
#endif
#endif
