#include "solver.h"
#include "utilities.h"
//////////////////////////////////////////////////////////////////////////
/* Computing terms needed to build the linear system solver   */
void compute_terms_for_matrix(cube_str *CUBE, int grid_points, double alpha, \
     double **inverse_term, double *inverse_term_level) {

     int i_atom, i_dim, i_grid;
     double delta_dist[NDIM], grid_pos_tmp[NDIM], atom_pos_tmp[NDIM]; 
     double Phi_all, Phi_recp, Phi_real;

     #pragma code_align 64
     #pragma omp parallel for private(i_atom,i_dim,i_grid,delta_dist, \
     grid_pos_tmp,atom_pos_tmp,Phi_real,Phi_recp,Phi_all)
     for(i_atom=0;i_atom<CUBE->n_atoms;i_atom++) {
        #pragma code_align 64
        #pragma unroll_and_jam
        for(i_dim=0;i_dim<NDIM;i_dim++) {
           atom_pos_tmp[i_dim] = CUBE->atom_pos[i_atom][i_dim];
        }
        #pragma code_align 64
        #pragma unroll_and_jam
        for(i_grid=0;i_grid<grid_points;i_grid++) {
           #pragma code_align 64
           #pragma unroll_and_jam
           for(i_dim=0;i_dim<NDIM;i_dim++) { 
              grid_pos_tmp[i_dim] = grid_pos[i_grid][i_dim];
              delta_dist[i_dim] = grid_pos_tmp[i_dim] - \
              atom_pos_tmp[i_dim];
           }
           if(fit_flag == 0) {
             Phi_all = Phi_coulomb(grid_pos_tmp,atom_pos_tmp);
           } else if(fit_flag == 1) {
             //Phi_recp = NEW_Phi_recp_coeff(delta_dist);
             Phi_recp = Phi_recp_coeff(delta_dist);
             //printf("Phi recp %e\n",Phi_recp);
             Phi_real = \
             Phi_real_coeff(alpha,grid_pos_tmp,atom_pos_tmp);
             //printf("Phi real %e\n",Phi_real);
             Phi_all = Phi_real + Phi_recp;
             //printf("Phi_all is %e\n",Phi_all);
           }
           inverse_term[i_grid][i_atom] = Phi_all;
           inverse_term_level[i_atom] += Phi_all;
        }
        inverse_term_level[i_atom] /= grid_points;
        //printf("term %d done\n",i_atom);       
     }    

}
///////////////////////////////////////////////////////////////////
/* Setting matrix and vector for linear solver */
void set_solver_data(cube_str *CUBE, int grid_points, double alpha) {

     int i_atom, j_atom, i_grid;
     double **inverse_term, *inverse_term_level, **matrix_trans, *vector_trans;

     // Memory allocation
     inverse_term = create_2d_double_array((long)grid_points,(long)CUBE->n_atoms,"For matrix");
     inverse_term_level = create_1d_double_array(CUBE->n_atoms,"For matrix");
     matrix_trans = create_2d_double_array((long)CUBE->n_atoms+1,(long)CUBE->n_atoms+1,"Temp matrix");
     vector_trans = create_1d_double_array(CUBE->n_atoms+1,"temp vector");     

     // Initializing temporary arrays 
     #pragma code_align 64
     #pragma unroll_and_jam
     #pragma ivdep
     for(j_atom=0;j_atom<CUBE->n_atoms;j_atom++) {
        vector_trans[j_atom] = 0.0;
        inverse_term_level[j_atom] = 0.0;
        #pragma code_align 64
        #pragma unroll_and_jam
        for(i_grid=0;i_grid<grid_points;i_grid++) {
           inverse_term[i_grid][j_atom] = 0.0;
        }
        #pragma code_align 64
        #pragma unroll_and_jam
        for(i_atom=0;i_atom<CUBE->n_atoms;i_atom++) {
           matrix_trans[i_atom][j_atom] = 0.0;
        }        
     }

     // Computing interaction terms
     compute_terms_for_matrix(CUBE,grid_points,alpha,inverse_term,\
     inverse_term_level);

     // Setting temporary Matrix and vector elements
     #pragma code_align 64
     #pragma unroll_and_jam
     for(i_atom=0;i_atom<CUBE->n_atoms;i_atom++) {
        #pragma code_align 64
        #pragma unroll_and_jam
        for(j_atom=i_atom;j_atom<CUBE->n_atoms;j_atom++) {
           #pragma code_align 64
           #pragma vector always
           for(i_grid=0;i_grid<grid_points;i_grid++) {
              // Matrix elements
              matrix_trans[i_atom][j_atom] += \
              (inverse_term[i_grid][j_atom] - inverse_term_level[j_atom])* \
              (inverse_term[i_grid][i_atom] - inverse_term_level[i_atom]);
              // Vector elements
              if(i_atom == j_atom) {
                vector_trans[i_atom] += V_pot_grid[i_grid]* \
                (inverse_term[i_grid][i_atom] - inverse_term_level[i_atom]);
              }
           }
           matrix_trans[j_atom][i_atom] = matrix_trans[i_atom][j_atom];
           //printf(" A nd B are %f %f\n",A_ik[i_atom][j_atom],B_k[i_atom]);
        }
     }
     // Setting the final linear system
     set_linear_matrix_vec(CUBE,matrix_trans,vector_trans);
     
     // Interfacing with CLAPACK
     integer lda, info;
     double *M;
     integer n=CUBE->n_atoms+1;
     integer ipiv[N_max];
     lda=1;

     M = create_1d_double_array(n*n,"");
     matrix_ctof(matrix_solv,n,n,M);
     // Solving the system
     dgesv_(&n,&lda,M,&n,ipiv,vector_solv,&n,&info);
     // Printing charge values
     if (info != 0) {
         printf("WARNING:  Linear solver reported an error:  eflag = %d\n\n",(int) info);
         printf("WARNING:  Charges may not erroneous ");
     }
     printf("\nFitted charges ordered as within the cube file\n");
     if(symm_flag == 1) {
	printf("   Symmetry applied:  only charges of symmetry unique charges printed");
     }
     printf("----------------------------------------------\n");


     for(i_atom=0;i_atom<CUBE->n_atoms;i_atom++) {
	printf("Charge %d of type %d = %f\n",i_atom+1,atom_index[i_atom],vector_solv[i_atom]);
     }

     printf("\n");

     if(symm_flag != 1) {
       const char* at_symbol[ 107 ]= { "Xx", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na",\
             "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe",\
             "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb",\
             "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba",\
             "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Th", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",\
             "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",\
             "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md",\
	     "No", "Lr", "Rf", "Db", "Sg" };

       printf("Atomic positions in Ang with Fitted charges \n");
       printf("  type        x (Ang)       y (Ang)       z (Ang)    charge\n");
       printf("------------------------------------------------------------\n");

       for(i_atom=0;i_atom<CUBE->n_atoms;i_atom++) {
  	printf("   %2s     %10.4f    %10.4f    %10.4f    %8.4f\n",at_symbol[CUBE->atom_index[i_atom]], \
                             CUBE->atom_pos[i_atom][0]*0.529177, CUBE->atom_pos[i_atom][1]*0.529177, 
  			   CUBE->atom_pos[i_atom][2]*0.529177, vector_solv[i_atom]);

       }

       printf("------------------------------------------------------------\n");
     }
     printf("\n"); 

     // Computing stats and printing QM and Coulomb ESPs on grid
     output_ESP_grid_data(inverse_term,inverse_term_level);

     // Releasing memory
     destroy_1d_double_array(M);
     destroy_1d_double_array(vector_trans);
     destroy_2d_double_array(matrix_trans);       
     destroy_1d_double_array(inverse_term_level);
     destroy_2d_double_array(inverse_term); 

}

