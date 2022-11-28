#include "vars_arrays.h"
#include "math.h"
#include "memory.h"
#include <complex>
#include "utilities.h"
/* Function prototypes */
double init_ewald_data(void);
void k_vectors_coeff(double);
double Phi_real_coeff(double,double *,double *);
double Phi_coulomb(double *,double *);
double Phi_recp_coeff(double *);
double NEW_Phi_recp_coeff(double *);
