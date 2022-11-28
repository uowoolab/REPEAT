#include "vars_arrays.h"
/////////////////////////////////////////////////////////
// Variables related to the input file
// Cube file name
char* cube_file_name;

// Fitting flag
int fit_flag;

// VDW scaling factors
double vdw_fact, vdw_fact_max;

// RESP fitting flag
int fit_RESP;

// Cutoff flag
int flag_cutoff;

// Cutoff radius
double R_cutoff;

// Symmmetry constrain flag
int symm_flag;

// QEq/Goddard flag restrain
int QEq_restraint;

// QEQ/Goddard restrain weight 
double lambda_prime;

// Total charge
double q_tot;

/////////////////////////////////////////////////////////
// Variables/arrays related to the Ewald sum calculation
int NMAX[NDIM];
double *KVEC;
int V_my_index_1[KMAX*KMAX+1], V_my_index_2[KMAX*KMAX+1];
int V_KSQ[KMAX+1][2*KMAX+1][2*KMAX+1];
// Real and reciprocal box vectors
double **real_box_vector, **recip_box_vector, Box_Volume;
////////////////////////////////////////////////////////
/* Arrays related to the VDW radii and QEq params */
double *vdw_radii, *elect_array, *chi_array;

///////////////////////////////////////////////////////
/* Arrays related to the RESP restrains */
double *str_array, *charge_eq_RESP;
///////////////////////////////////////////////////////
/* Arrays related to the grid ESP */
int counter_grid;
double **grid_pos, *V_pot_grid;
//////////////////////////////////////////////////////
/* Arrays related to the linear solver */
int *atom_index, n_atoms_all;
double **matrix_solv, *vector_solv;
//////////////////////////////////////////////////////
/* Arrays related to the symmetry restrains */
int **connectivity_array, *linking_array, *charges_idem, \
    *charges_equiv, order_type;

