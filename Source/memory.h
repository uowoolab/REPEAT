#include <stdlib.h>
#include <stdio.h>
/* Function prototypes */
void *smalloc(long int,const char *);
void sfree(void *);
double *create_1d_double_array(int,const char *); 
void destroy_1d_double_array(double *);
double **create_2d_double_array(long,long,const char *);
void destroy_2d_double_array(double **);
int *create_1d_int_array(int,const char *);
void destroy_1d_int_array(int *);
char *create_1d_char_array(int,const char *);
void destroy_1d_char_array(char *);
char **create_2d_char_array(int,int,const char *);
void destroy_2d_char_array(char **);
int **create_2d_int_array(int,int,const char *);
void destroy_2d_int_array(int **);

