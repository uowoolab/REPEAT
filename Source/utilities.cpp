#include "utilities.h"
using namespace std;
///////////////////////////////////////////////////////////////////
/* Converting from 2d to 1d indices */
int convert_2d_to_1d(int i, int j, int max) {

    int index;

    index = i + j*max;

    return index;
}
///////////////////////////////////////////////////////////////////
/* Relabeling indices in the reciprocal sum arrays  */
int relabel_index(int order) {

    int index;
  
    index = order + KMAX;

    return index;

}
//////////////////////////////////////////////////////////////////
/* Computing cross product of two vectors */
void cross_prod(double *vec_in_1, double *vec_in_2, double *vec_out) {

     vec_out[0] = vec_in_1[1]*vec_in_2[2] - vec_in_1[2]*vec_in_2[1];
     vec_out[1] = vec_in_1[2]*vec_in_2[0] - vec_in_1[0]*vec_in_2[2];
     vec_out[2] = vec_in_1[0]*vec_in_2[1] - vec_in_1[1]*vec_in_2[0];

}
/////////////////////////////////////////////////////////////////
/* Computing the dot product of two vectors */
double dot_prod(double *vec_in_1, double *vec_in_2) {

       int i_dim;
       double out;
 
       out = 0.0;  
       for(i_dim=0;i_dim<NDIM;i_dim++) {
          out += vec_in_1[i_dim]*vec_in_2[i_dim];
       }
       return out;   

}
//////////////////////////////////////////////////////////////////
/* Computing the real and reciprocal box vectors */
void compute_box_vectors(cube_str *CUBE) {

     int i_dim, j_dim;
     double *a, *b, *c, *vec_out;

     // Memory allocation
     real_box_vector = \
     create_2d_double_array((long)NDIM,(long)NDIM,"Real space box vectors");
     recip_box_vector = \
     create_2d_double_array((long)NDIM,(long)NDIM,"Reciprocal space box vectors");
     a = create_1d_double_array(NDIM,"tmp");
     b = create_1d_double_array(NDIM,"tmp");
     c = create_1d_double_array(NDIM,"tmp");
     vec_out = create_1d_double_array(NDIM,"tmp");

     // Building the real box vectors
     for(j_dim=0;j_dim<NDIM;j_dim++) {
        for(i_dim=0;i_dim<NDIM;i_dim++) {
           real_box_vector[j_dim][i_dim] = \
           CUBE->n_grid[j_dim]*CUBE->axis_vector[j_dim][i_dim];
        }
     }

     // Defining some temporal vectors
     a[0] = real_box_vector[0][0];
     a[1] = real_box_vector[0][1];
     a[2] = real_box_vector[0][2];
     b[0] = real_box_vector[1][0];
     b[1] = real_box_vector[1][1];
     b[2] = real_box_vector[1][2];
     c[0] = real_box_vector[2][0];
     c[1] = real_box_vector[2][1];
     c[2] = real_box_vector[2][2];

     // Computing the volume of the simulation box
     cross_prod(b,c,vec_out);
     Box_Volume = dot_prod(a,vec_out);  
     printf("Real box volume = %f bohrs^3 \n",Box_Volume);

     // Building the reciprocal box vectors
     cross_prod(a,b,vec_out);
     for(i_dim=0;i_dim<NDIM;i_dim++) {
        recip_box_vector[2][i_dim] = (2*PI/Box_Volume)*vec_out[i_dim];
     }
     cross_prod(c,a,vec_out);
     for(i_dim=0;i_dim<NDIM;i_dim++) {
        recip_box_vector[1][i_dim] = (2*PI/Box_Volume)*vec_out[i_dim];
     }
     cross_prod(b,c,vec_out);
     for(i_dim=0;i_dim<NDIM;i_dim++) {
        recip_box_vector[0][i_dim] = (2*PI/Box_Volume)*vec_out[i_dim];
     }

     // Releasing memory
     destroy_1d_double_array(a);
     destroy_1d_double_array(b);
     destroy_1d_double_array(c);
     destroy_1d_double_array(vec_out);

}
////////////////////////////////////////////////////////////////////////
/* Repositioning atoms within the box */
void reposition_atoms(cube_str *CUBE, double *a, double *b, double *c) {

     int i_atom, i_dim;
     double *r;

     // Memory allocation
     r = create_1d_double_array(NDIM,"tmp");

     // Setting origin of the real box to (0,0,0)
     for(i_atom=0;i_atom<CUBE->n_atoms;i_atom++) {
        for(i_dim=0;i_dim<NDIM;i_dim++) {
           CUBE->atom_pos[i_atom][i_dim] -= CUBE->axis_zero[i_dim];
           r[i_dim] = CUBE->atom_pos[i_atom][i_dim];
        }
        // Repositioning atoms inside the box if necessary (PBCs active)
        if(fit_flag == 1) {
           // First vector
           iterate(CUBE,i_atom,r,b,c,a);                      
           // Second vector
           iterate(CUBE,i_atom,r,c,a,b);
           // Third vector
           iterate(CUBE,i_atom,r,a,b,c);           
        }
     }

     // Releasing memory
     destroy_1d_double_array(r);

}
//////////////////////////////////////////////////////////////
/* Shifting atomic positions into box along vec_3 direction  */
void iterate(cube_str *CUBE, int i_atom, double *vec_r, double *vec_1, \
     double *vec_2, double *vec_3) {

     int i_dim;
     double *vec_out, dot, dot_e, proj;

     // Memory allocation
     vec_out = create_1d_double_array(NDIM,"tmp");

     cross_prod(vec_r,vec_1,vec_out);
     dot = dot_prod(vec_2,vec_out);
     cross_prod(vec_3,vec_1,vec_out);
     dot_e = dot_prod(vec_2,vec_out);
     proj = dot/dot_e;
     // Shift atomic positions according to the vectorial
     // projection computed above
     if(proj < 0.0) {
       for(i_dim=0;i_dim<NDIM;i_dim++) {
          CUBE->atom_pos[i_atom][i_dim] += vec_3[i_dim];
       }
     } else if(proj > 1.0) {
       for(i_dim=0;i_dim<NDIM;i_dim++) {
           CUBE->atom_pos[i_atom][i_dim] -= vec_3[i_dim];
       }
     }

     // Releasing memory
     destroy_1d_double_array(vec_out);

}
////////////////////////////////////////////////////////////////
/* Convert from C to fortran matrix form */
void matrix_ctof(double **in, int rows, int cols, double *out) {

   int i, j;

     for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i+j*cols] = in[i][j];

}
////////////////////////////////////////////////////////////////

void error_message(char *err) {

    printf("\n\n------*ERROR* REPEAT did not complete normally *Error*------\n");
    printf(err);
    printf("\n\n\n");
    exit(0);

}
