// Carlos Campana: Nov 2013
// REPEAT C++/C-vanilla: a la LAMMPS style ;-)
//
//
//
#include "io_processing.h"
#include "setup.h"
#include "ewald.h"
#include "solver.h"
using namespace std;                            
int main(int narg, char **arg) {

    int grid_points;
    double alpha;

    printf("===================================================================\n");
    printf("REPEAT  Electrostatic potential fitted charges for periodic systems\n");
    printf("            Version 2.0 updated: Nov. 2022\n\n");
    printf("            Carlos Campañá and Tom Woo \n\n");
    printf(" Please reference:\n");
    printf("    Campañá et al. J. Chem. Theory Comput. 2009, 5, 10, 2866–2878\n");
    printf("		    https://doi.org/10.1021/ct9003405\n\n");
    printf("===================================================================\n\n");

    cube_str* CUBE;
    CUBE = new cube_str;
    // I/O processing of input parameter file
    read_input_file();

    // Reading ESP cube file
    read_cube_file(CUBE);

    // Reading connectivity file
    read_connectivity_file(CUBE);

    // Tabulate VDW radii, electronegativities and hardnesses
    setup_tabulated_arrays(); 

    // Reading RESP data file if necessary
    read_RESP_data_file();

    // Reading QEq data file if necessary
    read_QEq_data_file();

    // Computing the real and reciprocal box vectors 
    compute_box_vectors(CUBE);

    // Setting up the ESP grid data
    grid_points = setup_grid_potential(CUBE);

    // Computing reciprocal space parameters
    alpha = init_ewald_data();

    printf("\n======================================================= \n");
    printf("charges being fit...\n");
    printf("    -this can take minutes to hours depending on the \n");
    printf("     size of the system and number of grid points.\n");

    // Set and solve linear system
    set_solver_data(CUBE, grid_points, alpha); 

    printf("\nNormal Termination of REPEAT\n\n");

    // Clean up memory
    release_memory_global(CUBE);
    delete CUBE;

}

