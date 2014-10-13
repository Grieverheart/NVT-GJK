#ifndef OCTAHEDRON_H
#define OCTAHEDRON_H
#include "../include/typedefs.h"
#ifdef _MAIN_
//Definitions for main
#define RC 1.8
static const double iscrb_d2=1.0;
static const double cscrb_d2=3.0;
static const uint TriadV1=0,TriadV2=5; //Which linearly independent Vertices to use. 0,5 for Cube. 1,2 for Octahedron
#else
//Definitions for GJK
static const uint vert_neigh[][3]={{1,3,5},{0,2,6},{1,3,7},{0,2,4},{3,5,7},{0,4,6},{1,5,7},{2,4,6}};

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
	//	hill-climbing algorithm							//
	//////////////////////////////////////////////////////
	uint next=0,last=0,curr=0;
	double p=0.0;
	double max=dot_p(vertices,axis);
	while(1){
		for(uint i=0;i<3;i++){
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
			if(i==2) return curr;
		}
	}
}

static inline uint poly_min(double *vertices, double *axis){
	//////////////////////////////////////////////////////
	//  Find the maximum projection on the axis using a //
	//	hill-climbing algorithm							//
	//////////////////////////////////////////////////////
	uint next=0,last=0,curr=0;
	double p=0.0;
	double min=dot_p(vertices,axis);
	while(1){
		for(uint i=0;i<3;i++){
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
			if(i==2)return curr;
		}
	}
}
#endif
#endif
