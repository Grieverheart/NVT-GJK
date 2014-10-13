#include <stdio.h>
#include "../include/quaternion.h"

void quat_rot(quat_t a, REAL* b, uint nVertices){
	///////////////////////////////////////////
	// Multiply vector b by quaternion a	 //
	///////////////////////////////////////////
	REAL cross_temp[3];
	
	for(uint i = 0; i < nVertices; i++){
		cross_temp[0] = a.el[2] * b[3*i+2] - a.el[3] * b[3*i+1] + a.el[0] * b[3*i+0];
		cross_temp[1] = a.el[3] * b[3*i+0] - a.el[1] * b[3*i+2] + a.el[0] * b[3*i+1];
		cross_temp[2] = a.el[1] * b[3*i+1] - a.el[2] * b[3*i+0] + a.el[0] * b[3*i+2];
		
		b[3*i+0] += 2.0 * (a.el[2] * cross_temp[2] - a.el[3] * cross_temp[1]);
		b[3*i+1] += 2.0 * (a.el[3] * cross_temp[0] - a.el[1] * cross_temp[2]);
		b[3*i+2] += 2.0 * (a.el[1] * cross_temp[1] - a.el[2] * cross_temp[0]);
	}
}
