#include "vars_arrays.h"
#include "memory.h"
#include "math.h"
#include <string>
/* Function prototypes */
int convert_2d_to_1d(int,int,int);
int relabel_index(int);
void cross_prod(double *,double *,double *);
double dot_prod(double *, double *);
void compute_box_vectors(cube_str *);
void reposition_atoms(cube_str *, double *, double *, double *);
void iterate(cube_str *, int, double *, double *, double *, double *);
void matrix_ctof(double**,int,int,double*);
void error_message(char *);
