#include "io_processing.h"
#include <string>
using namespace std;
/////////////////////////////////////////////////////////////////
/* Reading the REPEAT input parameter file if it exists,
   otherwise creating a default one */
void read_input_file(void) {

     char *file_name, *line;
     FILE *default_param_file;

     // Memory allocation
     cube_file_name = create_1d_char_array(MAXLINE,"cube file name");
     line = create_1d_char_array(MAXLINE,"just a line");
     file_name = create_1d_char_array(MAXLINE,"input param file");
     file_name = (char *) "repeat.input";

     default_param_file = fopen(file_name,"r");
     if(default_param_file != NULL) {       
       fgets(line,MAXLINE,default_param_file);       
       fgets(line,MAXLINE,default_param_file);
       // Read name of cube file
       sscanf(line, "%s\n",cube_file_name);
       fgets(line,MAXLINE,default_param_file);
       fgets(line,MAXLINE,default_param_file);
       // Read fitting flag
       sscanf(line, "%d\n",&fit_flag);       
       fgets(line,MAXLINE,default_param_file);
       fgets(line,MAXLINE,default_param_file);
       // Read VDW scaling factor
       sscanf(line, "%lg\n",&vdw_fact);
       fgets(line,MAXLINE,default_param_file);
       fgets(line,MAXLINE,default_param_file);
       // Read RESP flag
       sscanf(line, "%d\n",&fit_RESP);
       fgets(line,MAXLINE,default_param_file);
       fgets(line,MAXLINE,default_param_file);
       // Read cutoff flag
       sscanf(line, "%d\n",&flag_cutoff);
       fgets(line,MAXLINE,default_param_file);
       fgets(line,MAXLINE,default_param_file);
       // Read cutoff radius
       sscanf(line, "%lg\n",&R_cutoff);
       fgets(line,MAXLINE,default_param_file);
       fgets(line,MAXLINE,default_param_file);
       // Read symmetry restrain flag
       sscanf(line, "%d\n",&symm_flag);
       fgets(line,MAXLINE,default_param_file);
       fgets(line,MAXLINE,default_param_file);
       // Read QEq restraint flag
       sscanf(line, "%d\n",&QEq_restraint);
       fgets(line,MAXLINE,default_param_file);
       fgets(line,MAXLINE,default_param_file);
       // Read QEq restraint weights
       sscanf(line, "%lg\n",&lambda_prime);
       fgets(line,MAXLINE,default_param_file);
       fgets(line,MAXLINE,default_param_file);
       // Read vdw_fact_max
       sscanf(line, "%lg\n",&vdw_fact_max);
       fgets(line,MAXLINE,default_param_file);
       fgets(line,MAXLINE,default_param_file);
       // Read total charge
       sscanf(line, "%lg\n",&q_tot);
       fclose(default_param_file);

       printf("Input Settings read from: %s \n", file_name );
       printf("=========================================\n");
       printf(" Cube file: %s\n", cube_file_name);
       printf(" van der Waals scaling factor: %.2f\n", vdw_fact);
       if (fit_RESP > 0) {
           printf("     RESP-like harmonic restraints applied to charges \n");
       } 
       if (QEq_restraint > 0) {
           if (fit_RESP > 0) {
               error_message((char *)"QEq and RESP restraints should not be applied at the same time.");
           }
           printf(" QEq restraints applied to charges \n");
           printf("  weighting for QEq restraints: %.4f\n", lambda_prime);
       } 
       if (symm_flag > 0) {
           printf(" symmetry constraints imposed\n");
       } 
       else {
           printf(" no symmetry imposed\n");
       } 

       if (fit_flag > 0) {
           printf(" Periodic system calculation\n");
       } 
       else {
           printf(" Molecular system calculation\n");
       } 

       printf(" radial cut_off (Bohr): %.2f\n", R_cutoff);
       printf(" maximum vdW factor: %.2f\n", vdw_fact_max);
       printf(" net charge of system (e): %.2f\n", q_tot );
       printf("\n\n");
  
     } else {
       // Error: Repeat input file could not be found, write a default input file and quit

       printf("Could not open repeat.input file");
       default_param_file = fopen(file_name,"w");
       fprintf(default_param_file,"Input ESP file name in cube format\n");
       fprintf(default_param_file," filename.cube\n");
       fprintf(default_param_file,"Fit molecular(0) or periodic(1:default) system?\n");
       fprintf(default_param_file," 1\n");
       fprintf(default_param_file,"van der Waals scaling factor (default = 1.0)\n");
       fprintf(default_param_file," 0.90000\n");
       fprintf(default_param_file,"Apply RESP-like harmonic restraints?, no(0:default), yes(1)\n");
       fprintf(default_param_file," 0\n");
       fprintf(default_param_file,"Read cutoff radius? no(0), yes(1:default)\n");
       fprintf(default_param_file," 1\n");
       fprintf(default_param_file,"If flag above=1 provide R_cutoff next (in Bohrs)\n");
       fprintf(default_param_file," 20.00000\n");
       fprintf(default_param_file,"Apply symmetry constraints? no(0:default), yes(1)\n");
       fprintf(default_param_file," 0\n");
       fprintf(default_param_file,"Use QEq-type restraints? no(0:default), yes(1)\n");
       fprintf(default_param_file," 0\n");
       fprintf(default_param_file,"If flag above=1 then provide weight next\n");
       fprintf(default_param_file," 0.00000\n");
       fprintf(default_param_file,"van der Waals max scaling factor\n");
       fprintf(default_param_file," 6.0\n");
       fprintf(default_param_file,"System total charge\n");
       fprintf(default_param_file," 0.00\n");
       printf("  repeat.input file with default input parameters created.\n\n");
       fclose(default_param_file);

       error_message((char *)"repeat.input file not found. One with default values created.");


     }
     // Releasing memory
     destroy_1d_char_array(line);

}
////////////////////////////////////////////////////////////////////
/* Reading cube file */
void read_cube_file(cube_str* CUBE) {

     int i_atom, i_grid, j_grid, i_grid_2d, k_grid;
     double junk;
     char* line;
     FILE* cube_file;

     // Memory allocation
     line = create_1d_char_array(MAXLINE,"just a line");
     CUBE->n_grid = create_1d_int_array(NDIM,"grid resolution");
     CUBE->axis_zero = create_1d_double_array(NDIM,"origin vector");
     CUBE->axis_vector = create_2d_double_array((long)NDIM,(long)NDIM,"axis vectors");
 
     cube_file = fopen(cube_file_name,"r");
     if(cube_file != NULL) {
       // Skipping two lines
       fgets(line,MAXLINE,cube_file);
       fgets(line,MAXLINE,cube_file);
       // Reading total # of atoms and origin vector
       fgets(line,MAXLINE,cube_file);
       sscanf(line, "%d %lg %lg %lg\n",&CUBE->n_atoms,&CUBE->axis_zero[0], \
       &CUBE->axis_zero[1],&CUBE->axis_zero[2]);
       // Reading grid info and axis vectors
       fgets(line,MAXLINE,cube_file);
       sscanf(line, "%d %lg %lg %lg\n",&CUBE->n_grid[0], \
       &CUBE->axis_vector[0][0], &CUBE->axis_vector[0][1], \
       &CUBE->axis_vector[0][2]);
       fgets(line,MAXLINE,cube_file);
       sscanf(line, "%d %lg %lg %lg\n",&CUBE->n_grid[1], \
       &CUBE->axis_vector[1][0], &CUBE->axis_vector[1][1], \
       &CUBE->axis_vector[1][2]);
       fgets(line,MAXLINE,cube_file);
       sscanf(line, "%d %lg %lg %lg\n",&CUBE->n_grid[2], \
       &CUBE->axis_vector[2][0],&CUBE->axis_vector[2][1], \
       &CUBE->axis_vector[2][2]);


       // More memory allocation
       CUBE->atom_index = \
       create_1d_int_array(CUBE->n_atoms,"atom types");
       CUBE->atom_pos =  \
       create_2d_double_array((long)CUBE->n_atoms,(long)NDIM,"atomic positions");      
       CUBE->V_pot = create_2d_double_array((long)CUBE->n_grid[0]*CUBE->n_grid[1], \
       (long)CUBE->n_grid[2],"ESP grid tabulation");

       // Reading atomic info
       for(i_atom=0;i_atom<CUBE->n_atoms;i_atom++){
          fgets(line,MAXLINE,cube_file);
          sscanf(line,"%d %lg %lg %lg %lg\n",&CUBE->atom_index[i_atom],&junk, \
          &CUBE->atom_pos[i_atom][0],&CUBE->atom_pos[i_atom][1], \
          &CUBE->atom_pos[i_atom][2]);  
       } 
       // Reading ESP grid data 
       for(i_grid=0;i_grid<CUBE->n_grid[0];i_grid++) {
          for(j_grid=0;j_grid<CUBE->n_grid[1];j_grid++) {  
             i_grid_2d = convert_2d_to_1d(i_grid,j_grid,CUBE->n_grid[0]);
             for(k_grid=0;k_grid<CUBE->n_grid[2];k_grid++) {
                fscanf(cube_file,"%lE",&CUBE->V_pot[i_grid_2d][k_grid]);
                //cout << CUBE->V_pot[i_grid_2d][k_grid] << endl;
             }
          }
       }
       fclose(cube_file);


       //print coordinates of system read from cube:

       printf("System Details read from cube file\n");
       printf("==========================================\n");
       printf("  Cell Vectors (Bohr)\n");       
       printf("  -----------------------------\n");       
       printf("  %8.4f  %8.4f  %8.4f \n", CUBE->axis_vector[0][0], CUBE->axis_vector[0][1],\
		                       CUBE->axis_vector[0][2]);
       printf("  %8.4f  %8.4f  %8.4f \n", CUBE->axis_vector[1][0], CUBE->axis_vector[1][1],\
		                       CUBE->axis_vector[1][2]);
       printf("  %8.4f  %8.4f  %8.4f \n", CUBE->axis_vector[2][0], CUBE->axis_vector[2][1],\
		                       CUBE->axis_vector[2][2]);

       const char* at_symbol[ 107 ]= { "Xx", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na",\
               "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe",\
               "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb",\
               "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba",\
               "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Th", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",\
               "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",\
               "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",\
               "Md", "No", "Lr", "Rf", "Db", "Sg" };

       printf("  \n");       
       printf("  Atom positions (Bohr)\n");       
       printf("  -----------------------------------------------\n");       
       for(i_atom=0;i_atom<CUBE->n_atoms;i_atom++){
          printf("   %2s     %10.4f    %10.4f    %10.4f  \n",at_symbol[CUBE->atom_index[i_atom]], \
			   CUBE->atom_pos[i_atom][0],CUBE->atom_pos[i_atom][1], CUBE->atom_pos[i_atom][2]);
       } 

      printf("\n");

      printf("==========================================\n");

     } else {

       printf("ESP Cube file: %s could not be open!\n", cube_file_name);

       error_message((char *)"Cannot read your ESP Cube file!");
     }

     // Releasing memory
     destroy_1d_char_array(line);

} 
////////////////////////////////////////////////////////////////////
/* Reading RESP penalty data file */
void read_RESP_data_file(void) {

     int line_counter, i_atom;
     char *name, *line;
     FILE* resp_file;

     // Memory allocation
     line = create_1d_char_array(MAXLINE,"just a line");
     name = create_1d_char_array(MAXLINE,"element name");
     str_array = create_1d_double_array(VDW_types,"RESP weight");
     charge_eq_RESP = create_1d_double_array(VDW_types,"RESP charge");

     if(fit_RESP == 1) {
       resp_file = fopen("RESP_parameters.input","r");
       if(resp_file != NULL) {
         printf("RESP-like restraint parameters read from: RESP_parameters.input\n");
         line_counter = 0;
         while(fgets(line,MAXLINE,resp_file) != NULL) {
              sscanf(line, "%s %lg %lg\n",name,&str_array[line_counter], \
              &charge_eq_RESP[line_counter]);
              line_counter++;
         }
         fclose(resp_file);
       } else {
        
         printf("Error: Cannot read your RESP parameter file: RESP_parameters.input\n");
         resp_file = fopen("RESP_parameters.input","w");
         fprintf(resp_file,"%s %f %f\n", "H"  ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "He" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Li" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Be" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "B"  ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "C"  ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "N"  ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "O"  ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "F"  ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Ne" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Na" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Mg" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Al" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Si" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "P"  ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "S"  ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Cl" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Ar" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "K"  ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Ca" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Sc" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Ti" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "V"  ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Cr" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Mn" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Fe" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Co" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Ni" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Cu" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Zn" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Ga" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Ge" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "As" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Se" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Br" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Kr" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Rb" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Sr" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Y"  ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Zr" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Nb" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Mo" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Tc" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Ru" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Rh" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Pd" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Ag" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Cd" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "In" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Sn" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Sb" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Te" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "I"  ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Xe" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Cs" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Ba" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "La" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Ce" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Pr" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Nd" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Pm" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Sm" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Eu" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Gd" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Th" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Dy" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Ho" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Er" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Tm" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Yb" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Lu" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Hf" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Ta" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "W"  ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Re" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Os" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Ir" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Pt" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Au" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Hg" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Tl" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Pb" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Bi" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Po" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "At" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Rn" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Fr" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Ra" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Ac" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Th" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Pa" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "U"  ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Np" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Pu" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Am" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Cm" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Bk" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Cf" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Es" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Fm" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Md" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "No" ,0.0 ,0.0);
         fprintf(resp_file,"%s %f %f\n", "Lw" ,0.0 ,0.0);
         printf("RESP parameter file RESP_parameters.input with empty values has been created.\n");
         fclose(resp_file);

         error_message((char *)"RESP_parameters.input file not found. One with empty values created.");

       }

     } else {
       for(i_atom=0;i_atom<VDW_types;i_atom++) {   
          str_array[i_atom] = 0.0;
          charge_eq_RESP[i_atom] = 0.0;
       }

     }

     // Releasing memory
     destroy_1d_char_array(line);
     destroy_1d_char_array(name);
     
}

////////////////////////////////////////////////////////////////////
/* Reading QEq parameter file */
void read_QEq_data_file(void) {
      // electronegativities and hardness taken from Rappe and Goddard's QEq method 
      // hardness extracted from the Materials Studio (MS)
      // Keep in mind that MS tabulates the hardness which is J_00/2
      // UNITS ARE HARTREES


     int line_counter, i_atom;
     char *name, *line;
     FILE* QEq_file;

     // Memory allocation
     line = create_1d_char_array(MAXLINE,"just a line");
     name = create_1d_char_array(MAXLINE,"element name");
     elect_array = create_1d_double_array(QEq_types,"Electronegativities");
     chi_array = create_1d_double_array(QEq_types,"Hardnesses");


     if(QEq_restraint == 1) {
       QEq_file = fopen("QEq_parameters.input","r");
       if(QEq_file != NULL) {
         printf("QEq restraint parameters read from: QEq_parameters.input\n");
         line_counter = 0;
         while(fgets(line,MAXLINE,QEq_file) != NULL) {
              sscanf(line, "%s %lg %lg\n",name,&elect_array[line_counter], &chi_array[line_counter]);
              line_counter++;
         }
         fclose(QEq_file);
       } else {
        
         printf("Error: Cannot read your QEq parameter file: QEq_parameters.input\n");
         printf("QEq parameter file QEq_parameters.input with default values has been created.\n");
         printf(" parameters for H, V, Cu, Zn, C, N, O, F, Cl, Br, I taken from \n");
         printf(" DOI:10.1021/jz401479k and the rest from Materials Studio\n");
         QEq_file = fopen("QEq_parameters.input","w");
         fprintf(QEq_file,"%s %f %f\n", "H " ,0.166401 ,0.255225);
         fprintf(QEq_file,"%s %f %f\n", "He" ,0.449802 ,0.899603);
         fprintf(QEq_file,"%s %f %f\n", "Li" ,0.110466 ,0.175364);
         fprintf(QEq_file,"%s %f %f\n", "Be" ,0.179222 ,0.326547);
         fprintf(QEq_file,"%s %f %f\n", "B " ,0.149493 ,0.310892);
         fprintf(QEq_file,"%s %f %f\n", "C " ,0.199590 ,0.215240);
         fprintf(QEq_file,"%s %f %f\n", "N " ,0.245778 ,0.243355);
         fprintf(QEq_file,"%s %f %f\n", "O " ,0.320234 ,0.314869);
         fprintf(QEq_file,"%s %f %f\n", "F " ,0.235784 ,0.409058);
         fprintf(QEq_file,"%s %f %f\n", "Ne" ,0.396149 ,0.792298);
         fprintf(QEq_file,"%s %f %f\n", "Na" ,0.104476 ,0.168749);
         fprintf(QEq_file,"%s %f %f\n", "Mg" ,0.145193 ,0.271424);
         fprintf(QEq_file,"%s %f %f\n", "Al" ,0.111752 ,0.216449);
         fprintf(QEq_file,"%s %f %f\n", "Si" ,0.153168 ,0.256284);
         fprintf(QEq_file,"%s %f %f\n", "P " ,0.200757 ,0.293988);
         fprintf(QEq_file,"%s %f %f\n", "S " ,0.254594 ,0.329707);
         fprintf(QEq_file,"%s %f %f\n", "Cl" ,0.213918 ,0.267278);
         fprintf(QEq_file,"%s %f %f\n", "Ar" ,0.289578 ,0.579156);
         fprintf(QEq_file,"%s %f %f\n", "K " ,0.0889681 ,0.141114);
         fprintf(QEq_file,"%s %f %f\n", "Ca" ,0.118734 ,0.211671);
         fprintf(QEq_file,"%s %f %f\n", "Sc" ,0.135014 ,0.15111);
         fprintf(QEq_file,"%s %f %f\n", "Ti" ,0.142033 ,0.166471);
         fprintf(QEq_file,"%s %f %f\n", "V " ,0.150415 ,0.154972);
         fprintf(QEq_file,"%s %f %f\n", "Cr" ,0.136447 ,0.223945);
         fprintf(QEq_file,"%s %f %f\n", "Mn" ,0.169778 ,0.182273);
         fprintf(QEq_file,"%s %f %f\n", "Fe" ,0.152763 ,0.229605);
         fprintf(QEq_file,"%s %f %f\n", "Co" ,0.156402 ,0.233353);
         fprintf(QEq_file,"%s %f %f\n", "Ni" ,0.161252 ,0.237542);
         fprintf(QEq_file,"%s %f %f\n", "Cu" ,0.199512 ,0.127447);
         fprintf(QEq_file,"%s %f %f\n", "Zn" ,0.136010 ,0.164010);
         fprintf(QEq_file,"%s %f %f\n", "Ga" ,0.110209 ,0.220564);
         fprintf(QEq_file,"%s %f %f\n", "Ge" ,0.148868 ,0.252683);
         fprintf(QEq_file,"%s %f %f\n", "As" ,0.190651 ,0.27995);
         fprintf(QEq_file,"%s %f %f\n", "Se" ,0.236219 ,0.303616);
         fprintf(QEq_file,"%s %f %f\n", "Br" ,0.209178 ,0.321925);
         fprintf(QEq_file,"%s %f %f\n", "Kr" ,0.257239 ,0.514479);
         fprintf(QEq_file,"%s %f %f\n", "Rb" ,0.0856607 ,0.135675);
         fprintf(QEq_file,"%s %f %f\n", "Sr" ,0.111127 ,0.179333);
         fprintf(QEq_file,"%s %f %f\n", "Y " ,0.144017 ,0.165515);
         fprintf(QEq_file,"%s %f %f\n", "Zr" ,0.141408 ,0.208143);
         fprintf(QEq_file,"%s %f %f\n", "Nb" ,0.141004 ,0.216375);
         fprintf(QEq_file,"%s %f %f\n", "Mo" ,0.144164 ,0.2335);
         fprintf(QEq_file,"%s %f %f\n", "Tc" ,0.153425 ,0.236293);
         fprintf(QEq_file,"%s %f %f\n", "Ru" ,0.154711 ,0.23225);
         fprintf(QEq_file,"%s %f %f\n", "Rh" ,0.158349 ,0.233132);
         fprintf(QEq_file,"%s %f %f\n", "Pd" ,0.166618 ,0.222475);
         fprintf(QEq_file,"%s %f %f\n", "Ag" ,0.163016 ,0.23034);
         fprintf(QEq_file,"%s %f %f\n", "Cd" ,0.184992 ,0.290828);
         fprintf(QEq_file,"%s %f %f\n", "In" ,0.110135 ,0.205057);
         fprintf(QEq_file,"%s %f %f\n", "Sn" ,0.146516 ,0.229605);
         fprintf(QEq_file,"%s %f %f\n", "Sb" ,0.180031 ,0.245627);
         fprintf(QEq_file,"%s %f %f\n", "Te" ,0.213729 ,0.25915);
         fprintf(QEq_file,"%s %f %f\n", "I " ,0.199586 ,0.210207);
         fprintf(QEq_file,"%s %f %f\n", "Xe" ,0.222696 ,0.445392);
         fprintf(QEq_file,"%s %f %f\n", "Cs" ,0.080222 ,0.125753);
         fprintf(QEq_file,"%s %f %f\n", "Ba" ,0.10341 ,0.176099);
         fprintf(QEq_file,"%s %f %f\n", "La" ,0.113553 ,0.163898);
         fprintf(QEq_file,"%s %f %f\n", "Ce" ,0.11392 ,0.161693);
         fprintf(QEq_file,"%s %f %f\n", "Pr" ,0.11392 ,0.161693);
         fprintf(QEq_file,"%s %f %f\n", "Nd" ,0.11392 ,0.161693);
         fprintf(QEq_file,"%s %f %f\n", "Pm" ,0.11392 ,0.161693);
         fprintf(QEq_file,"%s %f %f\n", "Sm" ,0.11392 ,0.161693);
         fprintf(QEq_file,"%s %f %f\n", "Eu" ,0.11392 ,0.161693);
         fprintf(QEq_file,"%s %f %f\n", "Gd" ,0.11392 ,0.161693);
         fprintf(QEq_file,"%s %f %f\n", "Th" ,0.11392 ,0.161693);
         fprintf(QEq_file,"%s %f %f\n", "Dy" ,0.11392 ,0.161693);
         fprintf(QEq_file,"%s %f %f\n", "Ho" ,0.11392 ,0.161693);
         fprintf(QEq_file,"%s %f %f\n", "Er" ,0.11392 ,0.161693);
         fprintf(QEq_file,"%s %f %f\n", "Tm" ,0.11392 ,0.161693);
         fprintf(QEq_file,"%s %f %f\n", "Yb" ,0.11392 ,0.161693);
         fprintf(QEq_file,"%s %f %f\n", "Lu" ,0.11392 ,0.161693);
         fprintf(QEq_file,"%s %f %f\n", "Hf" ,0.139644 ,0.279289);
         fprintf(QEq_file,"%s %f %f\n", "Ta" ,0.178083 ,0.256064);
         fprintf(QEq_file,"%s %f %f\n", "W " ,0.171726 ,0.283551);
         fprintf(QEq_file,"%s %f %f\n", "Re" ,0.147729 ,0.284433);
         fprintf(QEq_file,"%s %f %f\n", "Os" ,0.180068 ,0.279289);
         fprintf(QEq_file,"%s %f %f\n", "Ir" ,0.198442 ,0.279289);
         fprintf(QEq_file,"%s %f %f\n", "Pt" ,0.206306 ,0.25621);
         fprintf(QEq_file,"%s %f %f\n", "Au" ,0.179847 ,0.190063);
         fprintf(QEq_file,"%s %f %f\n", "Hg" ,0.230413 ,0.305747);
         fprintf(QEq_file,"%s %f %f\n", "Tl" ,0.117595 ,0.213141);
         fprintf(QEq_file,"%s %f %f\n", "Pb" ,0.143319 ,0.259444);
         fprintf(QEq_file,"%s %f %f\n", "Bi" ,0.17235 ,0.274879);
         fprintf(QEq_file,"%s %f %f\n", "Po" ,0.154711 ,0.309422);
         fprintf(QEq_file,"%s %f %f\n", "At" ,0.174555 ,0.349111);
         fprintf(QEq_file,"%s %f %f\n", "Rn" ,0.197339 ,0.394679);
         fprintf(QEq_file,"%s %f %f\n", "Fr" ,0.073497 ,0.146994);
         fprintf(QEq_file,"%s %f %f\n", "Ra" ,0.104476 ,0.178892);
         fprintf(QEq_file,"%s %f %f\n", "Ac" ,0.126782 ,0.180068);
         fprintf(QEq_file,"%s %f %f\n", "Th" ,0.126782 ,0.180068);
         fprintf(QEq_file,"%s %f %f\n", "Pa" ,0.126782 ,0.180068);
         fprintf(QEq_file,"%s %f %f\n", "U " ,0.126782 ,0.180068);
         fprintf(QEq_file,"%s %f %f\n", "Np" ,0.126782 ,0.180068);
         fprintf(QEq_file,"%s %f %f\n", "Pu" ,0.126782 ,0.180068);
         fprintf(QEq_file,"%s %f %f\n", "Am" ,0.126782 ,0.180068);
         fprintf(QEq_file,"%s %f %f\n", "Cm" ,0.126782 ,0.180068);
         fprintf(QEq_file,"%s %f %f\n", "Bk" ,0.126782 ,0.180068);
         fprintf(QEq_file,"%s %f %f\n", "Cf" ,0.126782 ,0.180068);
         fprintf(QEq_file,"%s %f %f\n", "Es" ,0.126782 ,0.180068);
         fprintf(QEq_file,"%s %f %f\n", "Fm" ,0.126782 ,0.180068);
         fprintf(QEq_file,"%s %f %f\n", "Md" ,0.126782 ,0.180068);
         fprintf(QEq_file,"%s %f %f\n", "No" ,0.126782 ,0.180068);
         fprintf(QEq_file,"%s %f %f\n", "Lr" ,0.126782 ,0.180068);
         fprintf(QEq_file,"%s %f %f\n", "Rf" ,0.126782 ,0.180068);
         fprintf(QEq_file,"%s %f %f\n", "Db" ,0.126782 ,0.180068);
         fprintf(QEq_file,"%s %f %f\n", "Sg" ,0.126782 ,0.180068);
         fprintf(QEq_file,"%s %f %f\n", "Bh" ,0.126782 ,0.180068);
         fclose(QEq_file);

         error_message((char *)"QEq_parameters.input file not found. One with default values created.");

       }
     } else {
       for(i_atom=0;i_atom<QEq_types;i_atom++) {   
          elect_array[i_atom] = 0.0;
          chi_array[i_atom] = 0.0;
       }

     }

     // Releasing memory
     destroy_1d_char_array(line);
     destroy_1d_char_array(name);

}



////////////////////////////////////////////////////////////////////////
/* Reading connectivity-symmetry file*/
void read_connectivity_file(cube_str *CUBE) {

     int i_type, index, i_index, i_atom;
     char *line;
     FILE *connectivity_file;

     // Memory allocation
     line = create_1d_char_array(MAXLINE,"just a line");
     linking_array = create_1d_int_array(CUBE->n_atoms,"linked charges");
     charges_idem = create_1d_int_array(CUBE->n_atoms,"idem charges");
     charges_equiv = create_1d_int_array(CUBE->n_atoms,"equivalent charges");

     for(i_atom=0;i_atom<CUBE->n_atoms;i_atom++) {
        linking_array[i_atom] = 0;
        charges_idem[i_atom] = i_atom;
        charges_equiv[i_atom] = i_atom;
     }

     // Reading symmetry info
     if(symm_flag == 1) {
       connectivity_file = fopen("symmetry.input","r");
       if(connectivity_file != NULL) {
         fgets(line,MAXLINE,connectivity_file);
         sscanf(line, "%d\n",&order_type);        
         connectivity_array = create_2d_int_array(order_type,N_symm_max,"");
         for(i_type=0;i_type<order_type;i_type++) {
            fgets(line,MAXLINE,connectivity_file);
            fgets(line,MAXLINE,connectivity_file);
            sscanf(line, "%d\n",&index);
            connectivity_array[i_type][0] = index;
            for(i_index=1;i_index<index+1;i_index++) {
               fscanf(connectivity_file,"%d\n", \
               &connectivity_array[i_type][i_index]);
            }              
            for(i_index=2;i_index<index+1;i_index++) {          
               linking_array[connectivity_array[i_type][i_index]-1] = 1;
               charges_idem[connectivity_array[i_type][i_index]-1] = 
               connectivity_array[i_type][1]-1;
            }       
         }
         for(i_atom=0;i_atom<CUBE->n_atoms;i_atom++) {
            if(linking_array[i_atom] == 0) {
              charges_equiv[i_atom] = index_equiv(i_atom,linking_array);
            } else if(linking_array[i_atom] == 1) {
              charges_equiv[i_atom] = \
              index_equiv(charges_idem[i_atom],linking_array);
            }
         }
         fclose(connectivity_file);

       } else {
         printf("Error: Cannot read your input file used for symmetry: symmetry.input\n");
         connectivity_file = fopen("symmetry.input","w");
         fprintf(connectivity_file,"%d\n",2);
         fprintf(connectivity_file,"%s\n","# first tree");
         fprintf(connectivity_file,"%d\n",2);
         fprintf(connectivity_file,"%d %d\n",1,2);
         fprintf(connectivity_file,"%s\n","# second tree");
         fprintf(connectivity_file,"%d\n",1);
         fprintf(connectivity_file,"%d\n",3);
         fclose(connectivity_file);
         printf("A sample connectivity file symmetry.input for H2O has been created.\n ");
         printf("Modify given your needs.  See manual for more details\n");

         error_message((char *)"Symmetry constraints require a symmetry.input file.");
       }
     }
     //JPS -topmost now, try removing this....
     //destroy_1d_char_array(line); 

}
/////////////////////////////////////////////////////////////////////////
/* Printing stats and Coulomb potential on grid points */
void output_ESP_grid_data(double **inverse_term,double *inverse_term_level) {
  
     int i_grid, i_atom, j_atom;
     double *V_coul, avrg_Coul_ESP, avrg_QM_ESP, V_pot_sq_sum, sq_error;
     char *file_name;
     FILE *Coul_ESP_file;

     file_name = (char *)"Coul_ESP.dat";
     V_coul = create_1d_double_array(counter_grid,"");

     avrg_QM_ESP = 0.0;
     avrg_Coul_ESP = 0.0;
     sq_error = 0.0;
     V_pot_sq_sum = 0.0;
     Coul_ESP_file = fopen(file_name,"w");
     if(Coul_ESP_file != NULL) {
       for(i_grid=0;i_grid<counter_grid;i_grid++) {
          V_coul[i_grid] = 0.0;
          for(i_atom=0;i_atom<n_atoms_all;i_atom++) {
             j_atom  = charges_equiv[i_atom];
             V_coul[i_grid] += (inverse_term[i_grid][i_atom] - \
             inverse_term_level[i_atom])*vector_solv[j_atom];
          }
          sq_error += pow((V_coul[i_grid] - V_pot_grid[i_grid]),2);
          V_pot_sq_sum += pow(V_pot_grid[i_grid],2);
          avrg_Coul_ESP += V_coul[i_grid]/counter_grid;
          avrg_QM_ESP += V_pot_grid[i_grid]/counter_grid;
          fprintf(Coul_ESP_file,"%f %f %f %f %f\n",grid_pos[i_grid][0], \
          grid_pos[i_grid][1],grid_pos[i_grid][2],V_coul[i_grid],\
          V_pot_grid[i_grid]);
       }
       printf("\n"); 
       printf("ESP Coulomb average = %le\n",avrg_Coul_ESP);
       printf("ESP QM average = %le\n",avrg_QM_ESP);
       printf("RMS Error of the potential = %le\n",sqrt(sq_error/counter_grid));
       fclose(Coul_ESP_file);
     } else {

       error_message((char *)"Cannot write Coulomb ESP grid data to file: Coul_ESP.dat");

     }
     destroy_1d_double_array(V_coul);

}

