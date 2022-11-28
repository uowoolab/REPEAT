#include "vars_arrays.h"
#include "memory.h"
#include "ewald.h"
#include "utilities.h"
#include "f2c.h"
#include "clapack.h"
#include "symmetry.h"
#include "io_processing.h"
/* Function prototypes */
void compute_terms_for_matrix(cube_str *,int,double,double **,double *);
void set_solver_data(cube_str *,int,double);
