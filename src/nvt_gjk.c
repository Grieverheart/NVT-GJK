#define _MAIN_
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include "../include/typedefs.h"
#include "../include/quaternion.h"
#include "../include/common.h"
#include "../include/gjk.h"
#include "../include/boop.h"
//Mersenne Twister
#define DSFMT_MEXP 19937
#define HAVE_SSE2
#include "../include/dSFMT.h"
dsfmt_t dsfmt;
//Mersenne Twister
#define EMPTY -1

////////////////////////////*Global Variables*////////////////////////////////
/**/uint part_cells[P],n[D]; 									  		  /**/
/**/uint **cell_neigh;                                                    /**/
/**/REAL box[D] = {0}, Rcm[D] = {0};   	 							      /**/
/**/tPart particle[P];			   	 									  /**/
/**/int cell_list[P],*cell=NULL;    									  /**/
/**/uint nV;		//Number of Vertices		       	  				  /**/
/**/REAL *vertices;//Holds original vertices	 						  /**/
/**/uint mv_moves = 0;                                                    /**/
/**/uint rot_moves = 0;                                                   /**/
//////////////////////////////////////////////////////////////////////////////
//extern// const REAL iscrb_d2;						  				      /**/
//extern// const REAL cscrb_d2;						  				      /**/
//extern// const uint TriadV1,TriadV2;         						  	  /**/
//////////////////////////////////////////////////////////////////////////////

static void Triad(tPart p, REAL* restrict aa);

// Fast and accurate cube-root
static inline double cbrt_f(double k){
	const unsigned int B1 = 715094163;
	double t = 0.0;
	unsigned int* pt = (unsigned int*) &t;
	unsigned int* px = (unsigned int*) &k;
	pt[1] = px[1] / 3 + B1;
	double x = t;
	for(int i = 0; i < 2; i++){
		double tri = x*x*x;
		x = x * ((tri + k + k) / (tri + tri + k));
	}
	return x;
}

//~ static inline void pushFront(uint *array, uint index){
    //~ uint temp = array[index];
    //~ for(uint i = index; i >= 1; i--){
        //~ array[i] = array[i - 1];
    //~ }
    //~ array[0] = temp;
//~ }

static inline void pushFront(uint *array, uint index){
    uint temp = array[index];
    memmove(array+1, array, index*sizeof(uint));
    array[0] = temp;
}

static void parse_coords(const char *file, const char *delim){
	////////////////////////////////////////////////////////////////////////////
	//  Parse a file and extract the positions of the centers of the spheres  //
	////////////////////////////////////////////////////////////////////////////
	char line[128];
	char *t1=NULL;
	int i=1,u=0;
	tPart temp_particle[P];
	REAL invBoxMatrix[D] = {0};
	
	
	FILE *fp;
	fp=fopen(file,"r");
	if(!fp){
		printf("Couldn't open file %s\n",file);
		return;
	}
	
	printf("Parsing Coordinates in %s...\n",file);
	
	while(fgets(line,sizeof(line),fp)!= NULL){
		u=0;
		REAL rotation[4];
		if(line[strlen(line)-1]=='\n')line[strlen(line)-1]='\0'; // Remove the niewline from the end of the string
		for(t1=strtok(line,delim);t1!=NULL;t1=strtok(NULL,delim)){
			if(i>2){
				if(u<3) temp_particle[i-3].center[u]=atof(t1);
				else rotation[u-3] = atof(t1);
			}
			else if(i == 2 && (u % 4) == 0){
				box[u / 4] = atof(t1);
			}
			u++;
		}
		if(i > 2){
			quat_t q_rotation;
			rotation[0] *= PI/180.0;
			q_rotation.el[0] = cos(0.5*rotation[0]);
			q_rotation.el[1] = rotation[1] * sin(0.5*rotation[0]);
			q_rotation.el[2] = rotation[2] * sin(0.5*rotation[0]);
			q_rotation.el[3] = rotation[3] * sin(0.5*rotation[0]);
			quat_rot(q_rotation, particle[i-3].vertices, nV);
		}
		i++;
	}
	
	invBoxMatrix[0] = 1.0 / box[0];
	invBoxMatrix[1] = 1.0 / box[1];
	invBoxMatrix[2] = 1.0 / box[2];
    
	for(uint j = 0; j < P; j++){
		for(uint k = 0; k < D; k++){
			particle[j].center[k] = invBoxMatrix[k] * temp_particle[j].center[k];
            particle[j].real_center[k] = particle[j].center[k];
		}
	}
    
	printf("Finished Parsing Coordinates\n");
	free(t1);
	fclose(fp);
}

static void parse_obj(const char *file, const char *delim){
	//////////////////////////////////////////////////////////////////
	//  Parse an .obj file and initialize the particle properties.  //
	//	This is primarilly made for cubes.							//
	//////////////////////////////////////////////////////////////////
	char line[128];
	char *t1=NULL;
	int v=0,u=0;
	
	FILE *fp;
	fp=fopen(file,"r");
	if(!fp){
		printf("Couldn't open file %s\n",file);
		return;
	}
	
	printf("Parsing Particle Properties in %s...\n",file);
	vertices=malloc(3*40*sizeof(*vertices));//Change 40 if more than 40 vertices
	while(fgets(line,sizeof(line),fp)!= NULL){
		if(line[0]=='#'){
            if(fgets(line,sizeof(line),fp) == NULL) printf("Error parsing obj file\n");
        }
		u=0;
		if(line[strlen(line)-1]=='\n')line[strlen(line)-1]='\0'; // Remove the newline from the end of the string
		if(line[0]=='v'){
			for(t1=strtok(line,delim);t1!=NULL;t1=strtok(NULL,delim)){
				if(u>0)vertices[D*v+u-1]=atof(t1);
				u++;
			}
			v++;
		}
	}
	nV=v;
	printf("Done Parsing Properties\nCounted %d vertices\n",nV);
	printf("Initializing Particle Properties...\n");
	vertices=realloc(vertices,D*nV*sizeof(*vertices));
	for(uint k=0;k<P;k++){
		particle[k].index=k;
		
		particle[k].vertices=malloc(D*nV*sizeof(*(particle[0].vertices)));
		
		
		memcpy(particle[k].vertices, vertices, D*nV*sizeof(*vertices));
		
	}
	printf("Done Initializing.\n");
	free(t1);
	fclose(fp);
}

static void save_config(uint run, uint step){
	//////////////////////////////////////////////////////////////////
	// Writing the data at each parse is handled by this function.  //
	// The static int i makes sure the files are numbered.			//
	//////////////////////////////////////////////////////////////////

								
	char *string="Data/pid%d-run.%02u.%07u"; 
	char buffer[64];	
	sprintf(buffer, string, getpid(), run, step);
	
	FILE *fp;
	fp=fopen(buffer,"w");
	fprintf(fp,"%d\n",P);
	for(int j = 0; j < D * D - 1; j++){
        if(j % 4 == 0) fprintf(fp,"%f\t", box[j / 4]);
        else fprintf(fp,"0.0\t");
    }
	fprintf(fp,"%f\n",box[D-1]);
	
	for(uint i = 0; i < P; i++){
	    REAL real_center[3];
		for(uint j = 0; j < D; j++){
		    real_center[j] = box[j] * particle[i].center[j] - Rcm[j];
		    if(real_center[j] > box[j]) real_center[j] -= box[j];
		    else if(real_center[j] < 0.0) real_center[j] += box[j];
		    fprintf(fp,"%f\t", real_center[j]);
        }
		/* Print Angle Axis format */
		REAL Rot[4];
		Triad(particle[i],Rot);
		for(uint j=0;j<D;j++) fprintf(fp,"%f\t",Rot[j]);
		fprintf(fp,"%f\n",Rot[3]);
	}
	fclose(fp);
}

static inline uint mtrandi(uint a){
	/////////////////////////////////////////
	// Returns a random int from 0 to a-1  //
	/////////////////////////////////////////
    return dsfmt_genrand_uint32(&dsfmt)%a;
}

static inline REAL mtrandf(REAL a,uint b){
	/////////////////////////////////////////////////////////////////////////
	// Returns a random REAL from -a to a if b is 1 and 0 to a if b is 0 //
	/////////////////////////////////////////////////////////////////////////
    return b?a*(2*dsfmt_genrand_open_open(&dsfmt)-1):a*dsfmt_genrand_open_open(&dsfmt);
}

static inline REAL sqr_d(REAL a){
	//////////////////////////////////////
	//  Returns the square of a REAL  //
	//////////////////////////////////////
	return a*a;
}

static void Triad(tPart p, REAL* restrict aa){
	////////////////////////////////////////////////////
	// Find the Axis-Angle Representation using the   //
	// Triad method of attitude estimation.			  //
	// Algorithm for Angle-Axis conversion can be	  //
	// found in "Rotation: ..." by Rebecca M. Brannon //
	////////////////////////////////////////////////////
	REAL a[3], b[3], c[3], A[3], B[3], C[3], R[9];
	REAL norm_a = 0.0, norm_b = 0.0, norm_c = 0.0;
	REAL norm_A = 0.0, norm_B = 0.0, norm_C = 0.0;
	
	/* Get Orthonormal Basis */
	
	b[0]=vertices[3*TriadV1+1]*vertices[3*TriadV2+2]-vertices[3*TriadV1+2]*vertices[3*TriadV2+1];
	b[1]=vertices[3*TriadV1+2]*vertices[3*TriadV2]-vertices[3*TriadV1]*vertices[3*TriadV2+2];
	b[2]=vertices[3*TriadV1]*vertices[3*TriadV2+1]-vertices[3*TriadV1+1]*vertices[3*TriadV2];
	
	c[0]=vertices[3*TriadV1+1]*b[2]-vertices[3*TriadV1+2]*b[1];
	c[1]=vertices[3*TriadV1+2]*b[0]-vertices[3*TriadV1]*b[2];
	c[2]=vertices[3*TriadV1]*b[1]-vertices[3*TriadV1+1]*b[0];
	
	B[0]=p.vertices[3*TriadV1+1]*p.vertices[3*TriadV2+2]-p.vertices[3*TriadV1+2]*p.vertices[3*TriadV2+1];
	B[1]=p.vertices[3*TriadV1+2]*p.vertices[3*TriadV2]-p.vertices[3*TriadV1]*p.vertices[3*TriadV2+2];
	B[2]=p.vertices[3*TriadV1]*p.vertices[3*TriadV2+1]-p.vertices[3*TriadV1+1]*p.vertices[3*TriadV2];
	
	C[0]=p.vertices[3*TriadV1+1]*B[2]-p.vertices[3*TriadV1+2]*B[1];
	C[1]=p.vertices[3*TriadV1+2]*B[0]-p.vertices[3*TriadV1]*B[2];
	C[2]=p.vertices[3*TriadV1]*B[1]-p.vertices[3*TriadV1+1]*B[0];
	
	for(uint i=0;i<3;i++){
		norm_a += sqr_d(vertices[3*TriadV1+i]);
		norm_b += sqr_d(b[i]);
		norm_c += sqr_d(c[i]);
		norm_A += sqr_d(p.vertices[3*TriadV1+i]);
		norm_B += sqr_d(B[i]);
		norm_C += sqr_d(C[i]);
	}
	norm_a = sqrt(norm_a);
	norm_b = sqrt(norm_b);
	norm_c = sqrt(norm_c);
	norm_A = sqrt(norm_A);
	norm_B = sqrt(norm_B);
	norm_C = sqrt(norm_C);
	
	for(uint i=0;i<3;i++){
		a[i]=vertices[3*TriadV1+i]/norm_a;
		b[i]/=norm_b;
		c[i]/=norm_c;
		A[i]=p.vertices[3*TriadV1+i]/norm_A;
		B[i]/=norm_B;
		C[i]/=norm_C;
	}
	
	/* Calculate the Rotation Matrix */
	
	for(uint i=0;i<3;i++){
		for(uint j=0;j<3;j++) R[3*i+j]=A[i]*a[j]+B[i]*b[j]+C[i]*c[j];
	}
	/* Convert Rotation Matrix to Angle-Axis Representation */
	REAL t=R[0]+R[4]+R[8];
	REAL cos_a=0.5*(t-1);
	REAL diff_cos=sqrt(1-sqr_d(cos_a));
	aa[0]=acos(cos_a);
	if(diff_cos>0.1e-10){
		REAL ss=1.0/(2.0*sin(aa[0]));
		aa[1]=ss*(R[7]-R[5]);
		aa[2]=ss*(R[2]-R[6]);
		aa[3]=ss*(R[3]-R[1]);
	}
	else if(cos_a>0.0){
		aa[1]=1.0;
		aa[2]=0.0;
		aa[3]=0.0;
	}
	else{
		uint column=0;
		REAL norm;
		for(uint i=0;i<3;i++){
			norm=0.0;
			uint flag=0;
			for(uint j=0;j<3;j++){
				if(R[3*j+i]==0.0)flag+=1;
				else norm+=sqr_d(R[3*j+i]);
			}
			if(flag<3){
				column=i;
				break;
			}
		}
		norm=sqrt(norm);
		
		aa[0]=PI;
		aa[1]=R[column]/norm;
		aa[2]=R[3+column]/norm;
		aa[3]=R[6+column]/norm;
		
	}
	aa[0]=aa[0]*180.0/PI;
}

static inline REAL min_dist(
	const REAL* restrict target,
	const REAL* restrict origin,
	const REAL* restrict box_measure,
	REAL* restrict min_vec)
{
	///////////////////////////////////////////////////////////////////
	//  Find the mimnimum distance by using the boundary conditions  //
	///////////////////////////////////////////////////////////////////
	/* Calculate the distance vector according to the boundary conditions */
	for(uint i = 0; i < D; i++){
		min_vec[i] = target[i] - origin[i];
		if(min_vec[i] > 0.5) min_vec[i] -= 1.0;
		else if(min_vec[i] < -0.5) min_vec[i] += 1.0;
	}
	/* Convert the distance vector to real coords and find its' square */
	min_vec[0] *= box_measure[0];
	min_vec[1] *= box_measure[1];
	min_vec[2] *= box_measure[2];
	
	return sqr_d(min_vec[0]) + sqr_d(min_vec[1]) + sqr_d(min_vec[2]);
}

static inline uint cell_index(const REAL* restrict position){
	/////////////////////////////////////////////////////////////////
	// Calculate the cell index given the particle position.       //
	// The number of cells and the cell size in each direction	   //
	// are needed. Since we're using reduced coords r_cell=1.0/n   //
	/////////////////////////////////////////////////////////////////
	uint temp = 
        (uint)(position[0]*n[0]) +
			n[0]*((uint)(position[1]*n[1]) +
				n[1]*(uint)(position[2]*n[2]));
        
    return temp;
}

static inline uint neigh(uint icell, uint m){
	////////////////////////////////////////////////////////////////////////////////
	//  Give a neighbour of cell "icell" given the index "m". The number of cells //
	//	"n" for each dimension and the total number of cells is needed.			  //
	//	Note1: Cell icell (i.e. the current cell) occurs for m=4,13 for D=2,3	  //
	////////////////////////////////////////////////////////////////////////////////
	uint cell_d[D],new_cell_d[D];
	cell_d[0] = icell % n[0];
	cell_d[1] = (icell / n[0]) % n[1];
    cell_d[2] = icell / (n[0] * n[1]);
	
	new_cell_d[0] = (n[0] + cell_d[0] + m % 3 - 1) % n[0];
	new_cell_d[1] = (n[1] + cell_d[1] + (m / 3) % 3 - 1) % n[1];
    new_cell_d[2] = (n[2] + cell_d[2] + m / 9 - 1) % n[2];
	
    return new_cell_d[0] + (new_cell_d[1] + new_cell_d[2] * n[1]) * n[0];
}

static void create_list(void){ 
	////////////////////////////////////////////////////////////////////////
	//  Create a Cell linked list for the particles at "centers" in a box //
	//	with matrix "box". "cell" contains the cells and cell_list the 	  //
	//	particles in each cell.											  //
	////////////////////////////////////////////////////////////////////////
	
	uint ncells = 1;
	REAL vol = 1.0;
	
	for(uint i = 0; i < D; i++) vol *= box[i];
	
	for(uint i = 0; i < D; i++){
        n[i] = (uint)(box[i] / RC);
		ncells *= n[i];
	}
	
	cell = realloc(cell,ncells*sizeof(*cell));
	if(cell == NULL){
		printf("ERROR\n");
		save_config(0,666);
	}
	
    cell_neigh = malloc(ncells * sizeof(uint*));
	for(uint i = 0; i < ncells; i++){
        cell[i] = EMPTY; //Empty the cells
        cell_neigh[i] = malloc(27 * sizeof(uint));
        for(uint j = 0; j < 27; j++){
            cell_neigh[i][j] = neigh(i, j);
        }
    }
	uint icell;
	for(uint i = 0; i < P;i++){
		icell = cell_index(particle[i].center); //Calculate the Cell index. Note: Some pointer tricks ;) (Start from centers[D*i])
		cell_list[i] = cell[icell]; //Insert the particle in the lists
		cell[icell] = i;
		part_cells[i] = icell;
	}
}

static __attribute__ ((noinline)) void rotateParticle(REAL dtheta) { 
	///////////////////////////////////////////////////////////
	// This function attempts to rotate a random particle    //
	// about an axis passing through its' center.			 //
	///////////////////////////////////////////////////////////
	uint ra_sphere = mtrandi(P); //Pick a random particle
    const uint nb = 27; //Number of neighbours
	uint icell_ra = part_cells[ra_sphere];
	quat_t q_rotation;
	REAL vertices_old[D*nV];
	
	/* Sample Random Quaternion */
	
	/* Generate Random Unit Vector */
	/* Ref:Sphere Point Picking, Wolfram */
	REAL x1 = mtrandf(1.0,1);
	REAL x2 = mtrandf(1.0,1);
	REAL s1 = sqr_d(x1) + sqr_d(x2);
	while(s1 > 1.0){
		x1 = mtrandf(1.0,1);
		x2 = mtrandf(1.0,1);
		s1 = sqr_d(x1) + sqr_d(x2);
	}
	REAL s2 = sqrt(1.0 - s1);
	
	/* Generate Random Angle */
	REAL angle = mtrandf(dtheta * dtheta * dtheta, 0);
	angle = cbrt_f(angle);
	while(mtrandf(1.0,0) > sqr_d(sin(angle) / angle)){
		angle = mtrandf(dtheta * dtheta * dtheta, 0);
		angle = cbrt_f(angle);
	}
	REAL sina = sin(angle * 0.5);
	/* Assemble Random Quaternion */
	q_rotation.el[0] = cos(0.5 * angle);
	q_rotation.el[1] = 2 * x1 * s2 * sina;
	q_rotation.el[2] = 2 * x2 * s2 * sina;
	q_rotation.el[3] = (1.0 - 2.0 * s1) * sina;
	//Note: Maybe we should check if we need to normalize
	
	
	/* Save the old Vertex Positions */
	memcpy(vertices_old, particle[ra_sphere].vertices, D*nV*sizeof(*vertices_old));
	
	/* Rotate Vertices. Eventually, make this geometry-specific instead of rotating all the vertices.  */
	quat_rot(q_rotation, particle[ra_sphere].vertices, nV);
	
	/* Check for overlaps */
	for(uint i = 0; i < nb; i++){ //Loop through the neighbouring cells
		for(int a = cell[cell_neigh[icell_ra][i]]; a != EMPTY; a = cell_list[a]){ //This loop traverses the linked list
			if(a != ra_sphere){
				REAL dist_vec[D];
				REAL dist = min_dist(particle[ra_sphere].center, particle[a].center, box, dist_vec);
				if(dist < iscrb_d2){
					memcpy(particle[ra_sphere].vertices, vertices_old, D*nV*sizeof(*vertices_old));
                    if(cell_neigh[icell_ra][i] != cell_neigh[icell_ra][0]){ //Reorder neighboring cells for coherence
                        pushFront(cell_neigh[icell_ra], i);
                    }
					return;
				}
				else if(dist < cscrb_d2){
					bool overlap = gjk_overlap(&particle[ra_sphere], &particle[a], dist_vec);
					if(overlap){
						memcpy(particle[ra_sphere].vertices, vertices_old, D*nV*sizeof(*vertices_old));
                        if(cell_neigh[icell_ra][i] != cell_neigh[icell_ra][0]){ //Reorder neighboring cells for coherence
                            pushFront(cell_neigh[icell_ra], i);
                        }
						return;
					}
				}
			}
		}
	}
	/* Change Accepted */
    rot_moves++;
}

static __attribute__ ((noinline)) void moveParticle(REAL d) { 
	///////////////////////////////////////////////////////////
	// This function attempts to move a random particle at a //
	// random distance [-dx,dx]  							 //
	///////////////////////////////////////////////////////////
	uint ra_sphere = mtrandi(P); //Pick a random particle
    const uint nb = 27; //Number of neighbours
	uint icell_new; //The cell of the random particle
	REAL dx[D];
	REAL center_old[D];
	
	
	for(uint k = 0; k < D; k++){
		dx[k] = mtrandf(d, 1); //Pick a random displacement
		center_old[k] = particle[ra_sphere].center[k];
        REAL temp = 1.0 + particle[ra_sphere].center[k] + dx[k];
		particle[ra_sphere].center[k] = temp - (uint)(temp);
	}
	icell_new = cell_index(particle[ra_sphere].center); // Calcualate the cell index for the new position
	
	for(uint i = 0; i < nb; i++){ //Loop through the neighbouring cells
		for(int a = cell[cell_neigh[icell_new][i]]; a != EMPTY; a = cell_list[a]){ //This loop traverses the linked list
			if(a != ra_sphere){
				REAL dist_vec[D];
				REAL dist = min_dist(particle[ra_sphere].center, particle[a].center, box, dist_vec);
				if(dist < iscrb_d2){//Overlapping. Reset Center coordinates
					for(uint k = 0; k < D; k++) particle[ra_sphere].center[k] = center_old[k];
                    if(cell_neigh[icell_new][i] != cell_neigh[icell_new][0]){ //Reorder neighboring cells for coherence
                        pushFront(cell_neigh[icell_new], i);
                    }
					return;
				}
				else if(dist < cscrb_d2){
					bool overlap = gjk_overlap(&particle[ra_sphere], &particle[a], dist_vec);
					if(overlap){//Overlapping. Reset Center coordinates
						for(uint k = 0; k < D; k++) particle[ra_sphere].center[k] = center_old[k];
                        if(cell_neigh[icell_new][i] != cell_neigh[icell_new][0]){ //Reorder neighboring cells for coherence
                            pushFront(cell_neigh[icell_new], i);
                        }
						return;
					}
				}//Do Overlaps Check
			}
		}
	}
	for(uint k = 0; k < D; k++){
		particle[ra_sphere].real_center[k] += dx[k];
		Rcm[k] += (box[k] * dx[k]) / P;
	}
	/* Change Accepted. Update the linked lists */
	uint icell_old = part_cells[ra_sphere]; // The cell the particle resides in before updating
	if(icell_old != icell_new){ //The particle has changed cells so the linked lists need to be updated
		if(cell[icell_old] == ra_sphere)cell[icell_old] = cell_list[ra_sphere];
		else{
			for(int a = cell[icell_old]; a != EMPTY; a = cell_list[a]){
				if(cell_list[a] == ra_sphere){
					cell_list[a] = cell_list[ra_sphere];
					break;
				}
			}
		}
		cell_list[ra_sphere] = cell[icell_new]; //Insert the particle into its' new cell
		cell[icell_new] = ra_sphere;
		part_cells[ra_sphere] = icell_new;
	}
	mv_moves++;
}


int main(int argc,char *argv[]){
	
	////////////////////////////////////////	
	////    Simulation Parameters   	////
	////////////////////////////////////////
	/**/uint mcsteps = 10001;		/**/
	/**/uint eq_steps = 15000;		    /**/
	/**/uint substeps = P;				/**/
	/**/REAL d = 0.02;                  /**/
    /**/REAL dtheta =  0.12;	        /**/
    /**/uint wait_time = 10000;			/**/
	////////////////////////////////////////	
	
	unsigned long seed;
	seed = time(NULL);
	dsfmt_init_gen_rand(&dsfmt, 0 /*seed + getpid()*/);
	
	if(argc < 3){
		printf("Did not provide enough input\nExiting...\n");
		return 0;
	}
	
	parse_obj(argv[2]," ");
	parse_coords(argv[1],"\t");
	create_list();
	init_gjk_cache();
	
	uint indices[10] = {0};
	uint a[10] = {0};
    
    char *string = "Data/cryst%d.dat";
	char buffer[64];
	sprintf(buffer, string, getpid());
    FILE *fp;
    fp = fopen(buffer,"w");
    
	for(uint i = 0; i < mcsteps; i++){
    	for(uint j = 0; j < substeps; j++){
		    if(mtrandf(1.0, 0) < 0.5) moveParticle(d);
		    else rotateParticle(dtheta);
        }
            
        if((i + 1) % 100 == 0){
            REAL mv_ratio = mv_moves / (50.0 * P);
            REAL rot_ratio = rot_moves / (50.0 * P);
            if(mv_ratio > 0.3) d *= 1.02;
            else if(mv_ratio < 0.25) d *= 0.98;
            if(rot_ratio > 0.3) dtheta *= 1.02;
            else if(rot_ratio < 0.25) dtheta *= 0.98;
       		if(dtheta > 2.0 * PI) dtheta -= 2.0 * PI;
            mv_moves = 0.0;
            rot_moves = 0.0;
		    if((i+1)%10000 == 0){
		        double cryst = crystallinity(box, particle, 1.99, 0.7);
		        fprintf(fp,"%u\t%f\n", i, cryst);
		        fflush(fp);
				save_config(i + 1, 666);
		    }
        }

        //~ if(i >= eq_steps){
            //~ uint u = i - eq_steps;
            //~ static uint run = 0;
            //~ if(u % wait_time == 0 && run < 10) run++;
            //~ for(uint k = 0; k < run; k++){
                //~ if(u - k * wait_time == a[k]){
                    //~ save_config(k + 1, a[k]);
                    //~ uint anew = a[k];
                    //~ while(anew == a[k]){
		                //~ indices[k]++;
		                //~ a[k] = floor(pow(10.0, 7.0 * (indices[k] + 1) / 100.0) + 0.5) - 1;
                    //~ }
                //~ }
            //~ }
        //~ }
    }
    
	fclose(fp);
    
	for(uint k = 0; k < P; k++){
		free(particle[k].vertices);
	}
	free(cell);
	free(vertices);
    
    uint ncells = n[0] * n[1] * n[2];
	for(uint i = 0; i < ncells; i++) free(cell_neigh[i]);
    free(cell_neigh);
	free_gjk_cache();
	printf("Done\n");
	return 0;
}
