#include "setup.h"
#include "utilities.h"
using namespace std;
//////////////////////////////////////////////////////////////////////////////
/* Tabulating the VDW radii */
void setup_tabulated_arrays(void) {

     // Memory allocation
     vdw_radii = create_1d_double_array(VDW_types,"VDW radii");

     // the VDW radii for the elements taken from UFF ordered by atomic number
     // UNITS ARE IN BOHR
      vdw_radii[ 0 ]= 2.72687 ;
      vdw_radii[ 1 ]= 2.23177 ;
      vdw_radii[ 2 ]= 2.31586 ;
      vdw_radii[ 3 ]= 2.59365 ;
      vdw_radii[ 4 ]= 3.85788 ;
      vdw_radii[ 5 ]= 3.63867 ;
      vdw_radii[ 6 ]= 3.4582 ;
      vdw_radii[ 7 ]= 3.30702 ;
      vdw_radii[ 8 ]= 3.17852 ;
      vdw_radii[ 9 ]= 3.06419 ;
      vdw_radii[ 10 ]= 2.81853 ;
      vdw_radii[ 11 ]= 2.85443 ;
      vdw_radii[ 12 ]= 4.25094 ;
      vdw_radii[ 13 ]= 4.05819 ;
      vdw_radii[ 14 ]= 3.91835 ;
      vdw_radii[ 15 ]= 3.81252 ;
      vdw_radii[ 16 ]= 3.72937 ;
      vdw_radii[ 17 ]= 3.65473 ;
      vdw_radii[ 18 ]= 3.60182 ;
      vdw_radii[ 19 ]= 3.21159 ;
      vdw_radii[ 20 ]= 3.11332 ;
      vdw_radii[ 21 ]= 2.99994 ;
      vdw_radii[ 22 ]= 2.97065 ;
      vdw_radii[ 23 ]= 2.85632 ;
      vdw_radii[ 24 ]= 2.79774 ;
      vdw_radii[ 25 ]= 2.75144 ;
      vdw_radii[ 26 ]= 2.71365 ;
      vdw_radii[ 27 ]= 2.67774 ;
      vdw_radii[ 28 ]= 3.3023 ;
      vdw_radii[ 29 ]= 2.61066 ;
      vdw_radii[ 30 ]= 4.14133 ;
      vdw_radii[ 31 ]= 4.04401 ;
      vdw_radii[ 32 ]= 3.99677 ;
      vdw_radii[ 33 ]= 3.97315 ;
      vdw_radii[ 34 ]= 3.95803 ;
      vdw_radii[ 35 ]= 3.91268 ;
      vdw_radii[ 36 ]= 3.88717 ;
      vdw_radii[ 37 ]= 3.44025 ;
      vdw_radii[ 38 ]= 3.16057 ;
      vdw_radii[ 39 ]= 2.95175 ;
      vdw_radii[ 40 ]= 2.99049 ;
      vdw_radii[ 41 ]= 2.88372 ;
      vdw_radii[ 42 ]= 2.8327 ;
      vdw_radii[ 43 ]= 2.79963 ;
      vdw_radii[ 44 ]= 2.7675 ;
      vdw_radii[ 45 ]= 2.73916 ;
      vdw_radii[ 46 ]= 2.97443 ;
      vdw_radii[ 47 ]= 2.69097 ;
      vdw_radii[ 48 ]= 4.21692 ;
      vdw_radii[ 49 ]= 4.14984 ;
      vdw_radii[ 50 ]= 4.17629 ;
      vdw_radii[ 51 ]= 4.22354 ;
      vdw_radii[ 52 ]= 4.25188 ;
      vdw_radii[ 53 ]= 4.16118 ;
      vdw_radii[ 54 ]= 4.26795 ;
      vdw_radii[ 55 ]= 3.49883 ;
      vdw_radii[ 56 ]= 3.32781 ;
      vdw_radii[ 57 ]= 3.35993 ;
      vdw_radii[ 58 ]= 3.40718 ;
      vdw_radii[ 59 ]= 3.37789 ;
      vdw_radii[ 60 ]= 3.35143 ;
      vdw_radii[ 61 ]= 3.32592 ;
      vdw_radii[ 62 ]= 3.30041 ;
      vdw_radii[ 63 ]= 3.1823 ;
      vdw_radii[ 64 ]= 3.26072 ;
      vdw_radii[ 65 ]= 3.23899 ;
      vdw_radii[ 66 ]= 3.22104 ;
      vdw_radii[ 67 ]= 3.20403 ;
      vdw_radii[ 68 ]= 3.18797 ;
      vdw_radii[ 69 ]= 3.17002 ;
      vdw_radii[ 70 ]= 3.4393 ;
      vdw_radii[ 71 ]= 2.96781 ;
      vdw_radii[ 72 ]= 2.99522 ;
      vdw_radii[ 73 ]= 2.89978 ;
      vdw_radii[ 74 ]= 2.79113 ;
      vdw_radii[ 75 ]= 2.94797 ;
      vdw_radii[ 76 ]= 2.68341 ;
      vdw_radii[ 77 ]= 2.60215 ;
      vdw_radii[ 78 ]= 3.11143 ;
      vdw_radii[ 79 ]= 2.55585 ;
      vdw_radii[ 80 ]= 4.10732 ;
      vdw_radii[ 81 ]= 4.06008 ;
      vdw_radii[ 82 ]= 4.12905 ;
      vdw_radii[ 83 ]= 4.44936 ;
      vdw_radii[ 84 ]= 4.4881 ;
      vdw_radii[ 85 ]= 4.50227 ;
      vdw_radii[ 86 ]= 4.62983 ;
      vdw_radii[ 87 ]= 3.47426 ;
      vdw_radii[ 88 ]= 3.28623 ;
      vdw_radii[ 89 ]= 3.20875 ;
      vdw_radii[ 90 ]= 3.23521 ;
      vdw_radii[ 91 ]= 3.20781 ;
      vdw_radii[ 92 ]= 3.23521 ;
      vdw_radii[ 93 ]= 3.23521 ;
      vdw_radii[ 94 ]= 3.19458 ;
      vdw_radii[ 95 ]= 3.14261 ;
      vdw_radii[ 96 ]= 3.1549 ;
      vdw_radii[ 97 ]= 3.13033 ;
      vdw_radii[ 98 ]= 3.1171 ;
      vdw_radii[ 99 ]= 3.10482 ;
      vdw_radii[ 100 ]= 3.09348 ;
      vdw_radii[ 101 ]= 3.06892 ;
      vdw_radii[ 102 ]= 3.05758 ;

}
////////////////////////////////////////////////////////////////////////
/* Setting the valid ESP grid points according to the VDW radii */
int setup_grid_potential(cube_str* CUBE) {

     int i_atom, i_grid, j_grid, k_grid, i_grid_2d, flag_status, \
         i_dim, k_neigh, j_neigh, i_neigh, flag2;
     double sum_V_pot, sum_V_pot_grid, dist, delta_dist[NDIM], grid_tmp[NDIM];
     double *vdwRMAX;
     char *file_name;
     FILE *ESP_file;   
 
     file_name = "grid_ESP.dat";
     // Memory allocation
     V_pot_grid = \
     create_1d_double_array(CUBE->n_grid[0]*CUBE->n_grid[1]*CUBE->n_grid[2],"");
     grid_pos = \
     create_2d_double_array(CUBE->n_grid[0]*CUBE->n_grid[1]*CUBE->n_grid[2],NDIM,"");
     vdwRMAX = create_1d_double_array(CUBE->n_atoms,"");
     
     for(i_atom=0;i_atom<CUBE->n_atoms;i_atom++) {
        vdwRMAX[i_atom] = vdw_fact_max*vdw_radii[CUBE->atom_index[i_atom]-1];
     }

     ESP_file = fopen(file_name,"w");
     if(ESP_file != NULL) {
       // Constructing the grid and ESP on it
       sum_V_pot = 0.0;
       sum_V_pot_grid = 0.0;
       counter_grid = -1;
       #pragma code_align 64
       #pragma omp parallel for private(i_grid,j_grid,k_grid,flag_status,\
       i_dim,i_atom,i_grid_2d,grid_tmp,delta_dist,dist,k_neigh,j_neigh,flag2,\
       i_neigh) reduction(+:sum_V_pot,sum_V_pot_grid)
       for(k_grid=0;k_grid<CUBE->n_grid[2];k_grid++) {
          #pragma code_align 64
          for(j_grid=0;j_grid<CUBE->n_grid[1];j_grid++) {
             #pragma code_align 64
             #pragma vector always
             for(i_grid=0;i_grid<CUBE->n_grid[0];i_grid++) {
                i_grid_2d = convert_2d_to_1d(i_grid,j_grid,CUBE->n_grid[0]);
                sum_V_pot += CUBE->V_pot[i_grid_2d][k_grid];
                flag_status = 1;
                flag2 = 0;
                 #pragma code_align 64
                 #pragma vector always
                for(i_dim=0;i_dim<NDIM;i_dim++){
                   grid_tmp[i_dim] = i_grid*CUBE->axis_vector[0][i_dim] + \
                   j_grid*CUBE->axis_vector[1][i_dim] + \
                   k_grid*CUBE->axis_vector[2][i_dim];                
                }
                 
                #pragma code_align 64
                #pragma vector always
                for(i_atom=0;i_atom<CUBE->n_atoms;i_atom++) {
                   #pragma code_align 64
                   #pragma unroll_and_jam
                   for(k_neigh=-1;k_neigh<2;k_neigh++){
                      #pragma code_align 64
                      #pragma unroll_and_jam
                      for(j_neigh=-1;j_neigh<2;j_neigh++){
                         #pragma code_align 64
                         #pragma unroll_and_jam
                         for(i_neigh=-1;i_neigh<2;i_neigh++){
                            dist = 0.0;                     
                            #pragma code_align 64
                            #pragma unroll_and_jam
                            for(i_dim=0;i_dim<NDIM;i_dim++){
                               delta_dist[i_dim] = grid_tmp[i_dim] - \
                               (CUBE->atom_pos[i_atom][i_dim] + \
                                i_neigh*real_box_vector[0][i_dim] + \
                                j_neigh*real_box_vector[1][i_dim] + \
                                k_neigh*real_box_vector[2][i_dim]);
                                dist += pow(delta_dist[i_dim],2);
                            }
                            dist = sqrt(dist);
                            if(dist <= vdw_fact* \
                              vdw_radii[CUBE->atom_index[i_atom]-1]) {
                              flag_status = 0;
                              goto next_grid_point;
                            }
                            if(dist <= vdwRMAX[i_atom]) {
                              flag2 = 1;
                            }
                         }
                      }
                   }
                }
                if(flag2 == 0) {
                  flag_status = 0;
                  goto next_grid_point;
                }
                // Storing ESP data corresponding to valid grid points
                if(flag_status == 1) {
                  #pragma omp critical 
                  {
                  counter_grid++;
                  #pragma code_align 64
                  #pragma vector always
                  for(i_dim=0;i_dim<NDIM;i_dim++){
                     grid_pos[counter_grid][i_dim] = grid_tmp[i_dim];
                  }
                  V_pot_grid[counter_grid] = CUBE->V_pot[i_grid_2d][k_grid];
                  sum_V_pot_grid += V_pot_grid[counter_grid];
                  }
                  //cout << sum_V_pot_grid << endl;
                }
                next_grid_point:;
                // Printing grid ESP data 
                fprintf(ESP_file,"%f %f %f %d %f\n",grid_tmp[0],grid_tmp[1], \
                grid_tmp[2],flag_status,CUBE->V_pot[i_grid_2d][k_grid]); 
             }
             fprintf(ESP_file,"%s\n","	");
          }
       }
       fclose(ESP_file);
       counter_grid++;
       printf("Total # of valid grid points = %d \n",counter_grid); 
       printf("ESP sum over full grid = %le \n",sum_V_pot);
       printf("ESP sum over valid grid points = %le \n",sum_V_pot_grid);

       // Shifting the zero level of the ESP
       #pragma code_align 64
       #pragma vector always
       for(i_grid=0;i_grid<counter_grid;i_grid++) {
          V_pot_grid[i_grid] -= sum_V_pot_grid/counter_grid;
       }

     } else {

       char* message = "Cannot write ESP grid data!!!";
       error_message(message);


     }
//JPS ----topmost routine now, trying to ski[ this
     //destroy_1d_double_array(vdwRMAX);
     return counter_grid;

}
//////////////////////////////////////////////////////////////////
/* Releasing all global memory used */
void release_memory_global(cube_str* CUBE) {

     destroy_1d_int_array(CUBE->n_grid);
     destroy_1d_double_array(CUBE->axis_zero);
     destroy_2d_double_array(CUBE->axis_vector);
     destroy_1d_int_array(CUBE->atom_index);
     destroy_2d_double_array(CUBE->atom_pos);
     destroy_2d_double_array(CUBE->V_pot);
     destroy_1d_char_array(cube_file_name);   
     destroy_1d_double_array(vector_solv);
     destroy_2d_double_array(matrix_solv);
     destroy_1d_double_array(KVEC); 
     destroy_2d_double_array(grid_pos);
     destroy_1d_double_array(V_pot_grid);
     destroy_2d_double_array(real_box_vector);
     destroy_2d_double_array(recip_box_vector);
     destroy_1d_double_array(charge_eq_RESP);
     destroy_1d_double_array(str_array);
     destroy_1d_double_array(vdw_radii);
     destroy_1d_double_array(elect_array);     
     destroy_1d_double_array(chi_array);
     destroy_1d_int_array(charges_equiv);
     destroy_1d_int_array(charges_idem);
     destroy_1d_int_array(linking_array);  
     destroy_1d_int_array(atom_index);
     if(symm_flag == 1) {
       destroy_2d_int_array(connectivity_array);
     }

}
