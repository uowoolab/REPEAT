#include <string>
#include "constants.h"
////////////////////////////////////////////////////////////////////////
// Cube structure definition
#ifndef VARS_ARRAYS_H
#define VARS_ARRAYS_H
struct cube_str
{
   double **atom_pos; // Atomic positions
   int *atom_index; // Atomic indices
   double **axis_vector; // Vector origin and axis vectors
   int n_atoms, *n_grid; // Total # of atoms in the system
                         // and grid resolution
   double **V_pot; // ESP tabulated data
   double *axis_zero;
};
#endif
//////////////////////////////////////////////////////////////////////
// Variables related to the input file
// Cube file name
extern char* cube_file_name;

// Fitting flag
extern int fit_flag;

// VDW scaling factors
extern double vdw_fact, vdw_fact_max;

// RESP fitting flag
extern int fit_RESP;

// Cutoff flag
extern int flag_cutoff;
 
// Cutoff radius
extern double R_cutoff;

// Symmmetry constrain flag
extern int symm_flag;

// QEq/Goddard flag restraint
extern int QEq_restraint;

// Goddard restraint weight 
extern double lambda_prime;

// Total charge
extern double q_tot;

//////////////////////////////////////////////////////////
// Variables/arrays related to the Ewald sum calculation
extern int NMAX[NDIM];
extern int V_my_index_1[KMAX*KMAX+1], V_my_index_2[KMAX*KMAX+1];
extern int V_KSQ[KMAX+1][2*KMAX+1][2*KMAX+1];
extern double *KVEC;
// Real and reciprocal box vectors
extern double **real_box_vector, **recip_box_vector, Box_Volume;
/////////////////////////////////////////////////////////
/* Arrays related to the VDW radii and QEq params */
extern double *vdw_radii, *elect_array, *chi_array;
///////////////////////////////////////////////////////
/* Arrays related to the RESP restrains */
extern double *str_array, *charge_eq_RESP;
//////////////////////////////////////////////////////
/* Arrays related to the grid ESP */
extern int counter_grid;
extern double **grid_pos, *V_pot_grid;
//////////////////////////////////////////////////////
/* Arrays related to the linear solver */
extern int *atom_index, n_atoms_all;
extern double **matrix_solv, *vector_solv;
/////////////////////////////////////////////////////
/* Arrays related to the symmetry restrains */
extern int **connectivity_array, *linking_array, *charges_idem, \
           *charges_equiv, order_type;

