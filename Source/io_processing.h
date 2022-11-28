#include "vars_arrays.h"
#include <iostream>
#include "memory.h"
#include "utilities.h"
#include <stdio.h>
#include "math.h"
#include "constants.h"
#include "symmetry.h"
/* Function prototypes */
void read_input_file(void);
void read_cube_file(cube_str*);
void read_RESP_data_file(void);
void read_QEq_data_file(void);
void read_connectivity_file(cube_str*);
void output_ESP_grid_data(double**,double*);
