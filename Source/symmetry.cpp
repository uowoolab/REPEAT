#include "symmetry.h"
/////////////////////////////////////////////////////////////////////////
/* Computing the equivalent atomic indices */
int index_equiv(int index, int *linker_array) {

     int i_index, counter_sum, index_out;

     counter_sum = 0;
     for(i_index=0;i_index<index;i_index++) {
        counter_sum += linker_array[i_index];
     }
     index_out = index - counter_sum;
     return index_out;

}
//////////////////////////////////////////////////////////////////////////
/* Setting the matrix and vector for the linear solver */
void set_linear_matrix_vec(cube_str* CUBE, double **matrix_trans, double *vector_trans) {

    int i_type, j_type, type_list, i_atom, j_atom, index, counter;

    double **matrix_tmp, *weights_tmp, *vector_tmp, *weights;

    // Memory allocation
    weights = create_1d_double_array(CUBE->n_atoms,"");
    weights_tmp = create_1d_double_array(CUBE->n_atoms,"");
    matrix_tmp = create_2d_double_array(CUBE->n_atoms+1,CUBE->n_atoms+1,"");
    vector_tmp = create_1d_double_array(CUBE->n_atoms+1,"");

    for(j_atom=0;j_atom<CUBE->n_atoms;j_atom++) {
        weights[j_atom] = 1.0;
    }
    n_atoms_all = CUBE->n_atoms;
 
    // Resizing matrix and vector if needed
    if(symm_flag == 1) {
      for(i_type=0;i_type<order_type;i_type++) {
         type_list = connectivity_array[i_type][0];
         i_atom = connectivity_array[i_type][1]-1;
         #pragma ivdep
         for(j_type=0;j_type<type_list-1;j_type++) {
            index = connectivity_array[i_type][j_type+2]-1;
            vector_trans[i_atom] += vector_trans[index];
            weights[i_atom] += 1.0;
            for(j_atom=0;j_atom<CUBE->n_atoms;j_atom++) {
               matrix_trans[j_atom][i_atom] += \
               matrix_trans[j_atom][index];
            }
         }
      }
      for(i_type=0;i_type<order_type;i_type++) {
         type_list = connectivity_array[i_type][0];
         i_atom = connectivity_array[i_type][1]-1;
         #pragma ivdep
         for(j_type=0;j_type<type_list-1;j_type++) {
            index = connectivity_array[i_type][j_type+2]-1;
            #pragma ivdep
            for(j_atom=0;j_atom<CUBE->n_atoms;j_atom++) {
               matrix_trans[i_atom][j_atom] += \
               matrix_trans[index][j_atom];
            }
         }
      }
      // Finding independent atomic sites 
      counter = 0;
      for(i_atom=0;i_atom<CUBE->n_atoms;i_atom++) {
         if(linking_array[i_atom] == 0) counter++;
      }     
      // Creating temporary vector and matrix
      atom_index = create_1d_int_array(counter+1,"atom types");       
      counter = 0;
      for(i_atom=0;i_atom<CUBE->n_atoms;i_atom++) {
         if(linking_array[i_atom] == 0) {
           weights_tmp[counter] = weights[i_atom];
           atom_index[counter] = CUBE->atom_index[i_atom];
           vector_tmp[counter] = vector_trans[i_atom];
           for(j_atom=0;j_atom<CUBE->n_atoms;j_atom++) {
              matrix_tmp[j_atom][counter] = matrix_trans[j_atom][i_atom];
           }
           counter++;
         }
      }
      counter = 0;
      for(i_atom=0;i_atom<CUBE->n_atoms;i_atom++) {
         if(linking_array[i_atom] == 0) {
           for(j_atom=0;j_atom<CUBE->n_atoms;j_atom++) {
              matrix_tmp[counter][j_atom] = matrix_tmp[i_atom][j_atom];
           }
           counter++;
         }
      }
      CUBE->n_atoms = counter;
    } else if(symm_flag == 0) {
      atom_index = create_1d_int_array(CUBE->n_atoms,"atom types"); 
      for(i_atom=0;i_atom<CUBE->n_atoms;i_atom++) {
         atom_index[i_atom] = CUBE->atom_index[i_atom];
      }
    }

    matrix_solv = create_2d_double_array(CUBE->n_atoms+1,CUBE->n_atoms+1,"Matrix for solver");
    vector_solv = create_1d_double_array(CUBE->n_atoms+1,"R.H.S vector for solver");
    // Loading the final matrix and vector for the linear solver
    #pragma ivdep
    for(i_atom=0;i_atom<CUBE->n_atoms;i_atom++) {
       if(symm_flag == 1) {
         vector_solv[i_atom] = vector_tmp[i_atom];
       } else if(symm_flag == 0) {
         vector_solv[i_atom] = vector_trans[i_atom];
       }
       #pragma ivdep
       for(j_atom=i_atom;j_atom<CUBE->n_atoms;j_atom++) {
          if(symm_flag == 1) {
             matrix_solv[i_atom][j_atom] = matrix_tmp[i_atom][j_atom];
             matrix_solv[j_atom][i_atom] = matrix_tmp[j_atom][i_atom];
          } else if(symm_flag == 0) {
             matrix_solv[i_atom][j_atom] = matrix_trans[i_atom][j_atom];
             matrix_solv[j_atom][i_atom] = matrix_trans[j_atom][i_atom];
          }
       }
    }
    #pragma ivdep
    for(i_atom=0;i_atom<CUBE->n_atoms;i_atom++) {
       // RESP constraints
       matrix_solv[i_atom][i_atom] += \
       fit_RESP*str_array[atom_index[i_atom]-1];       
       if(symm_flag == 1) {
          matrix_solv[i_atom][CUBE->n_atoms] = weights_tmp[i_atom];
          if(QEq_restraint == 1) {
            matrix_solv[i_atom][i_atom] += \
            0.5*lambda_prime*weights_tmp[i_atom]* \
            chi_array[atom_index[i_atom]-1];
            vector_solv[i_atom] -= 0.5*lambda_prime*weights_tmp[i_atom]* \
            elect_array[atom_index[i_atom]-1];      
          }
       } else if(symm_flag == 0) {
          matrix_solv[i_atom][CUBE->n_atoms] = weights[i_atom]; 
          if(QEq_restraint == 1) {
            matrix_solv[i_atom][i_atom] += \
            0.5*lambda_prime*weights[i_atom]* \
            chi_array[atom_index[i_atom]-1];
            vector_solv[i_atom] -= 0.5*lambda_prime*weights[i_atom]* \
            elect_array[atom_index[i_atom]-1];
          }
       }
       matrix_solv[CUBE->n_atoms][i_atom] = \
       matrix_solv[i_atom][CUBE->n_atoms];
       vector_solv[i_atom] += \
       fit_RESP*str_array[atom_index[i_atom]-1]* \
       charge_eq_RESP[atom_index[i_atom]-1];
    }
    matrix_solv[CUBE->n_atoms][CUBE->n_atoms] = 0.0;
    vector_solv[CUBE->n_atoms] = q_tot;
  
    // Releasing memory
    destroy_1d_double_array(weights);
    destroy_1d_double_array(weights_tmp);
    destroy_1d_double_array(vector_tmp);
    destroy_2d_double_array(matrix_tmp);

}
