/* Memory allocation */
#include "memory.h"
////////////////////////////////////////////////////////////////////////
void *smalloc(long int order, const char *name) {

     if (order == 0) return NULL;
       void *ptr = malloc(order);
     if (ptr == NULL) {
        char str[128];
        sprintf(str,"Failed to allocate %ld bytes for array %s",order,name);
     }
     return ptr;

}
////////////////////////////////////////////////////////////////////////
void sfree(void *ptr) {

     if (ptr == NULL) return;
       free(ptr);

}
////////////////////////////////////////////////////////////////////////
/* Dynamic allocation of a 1d double array */
double *create_1d_double_array(int order, const char *name) {

       double *array = (double *) smalloc(order*sizeof(double),name);
       return array;

}
////////////////////////////////////////////////////////////////////////
/* Deallocation of a 1d double array */
void destroy_1d_double_array(double *array) {

     if (array == NULL) {
        printf("In destroy_1d_double_array pointer is NULL\n");
        return;
     }
     sfree(array);

}
////////////////////////////////////////////////////////////////////////
/* Dynamic allocation of a 2d double array */
double **create_2d_double_array(int n1, int n2, const char *name) {
       
       double *data = (double *) smalloc(n1*n2*sizeof(double),name);
       double **array = (double **) smalloc(n1*sizeof(double *),name);

       int order = 0, i;
       #pragma code_align 64
       for (i = 0; i < n1; i++) {
           array[i] = &data[order];
           order += n2;
       }
       return array;

}
////////////////////////////////////////////////////////////////////////
/* Deallocation of a 2d double array */
void destroy_2d_double_array(double **array) {

     if (array == NULL) {
        printf("In destroy_2d_double_array pointer is NULL\n");
        return;
     }
     sfree(array[0]);
     sfree(array);

}
////////////////////////////////////////////////////////////////////////
/* Dynamic allocation of a 1d integer array */
int *create_1d_int_array(int order, const char *name) {

       int *array = (int *) smalloc(order*sizeof(int),name);
       return array;

}
///////////////////////////////////////////////////////////////////////
/* Dynamic allocation of a 2d integer array */
int **create_2d_int_array(int n1, int n2, const char *name) {

       int *data = (int *) smalloc(n1*n2*sizeof(int),name);
       int **array = (int **) smalloc(n1*sizeof(int *),name);

       int order = 0, i;
       for (i = 0; i < n1; i++) {
           array[i] = &data[order];
           order += n2;
       }
       return array;

}
////////////////////////////////////////////////////////////////////////
/* Deallocation of a 2d double array */
void destroy_2d_int_array(int **array) {

     if (array == NULL) {
        printf("In destroy_2d_int_array pointer is NULL\n");
        return;
     }
     sfree(array[0]);
     sfree(array);

}
////////////////////////////////////////////////////////////////////////
/* Deallocation of a 1d integer array */
void destroy_1d_int_array(int *array) {

     if (array == NULL) {
        printf("In destroy_1d_int_array pointer is NULL\n");
        return;
     }
     sfree(array);

}
////////////////////////////////////////////////////////////////////////
/* Dynamic allocation of a 1d character array */
char *create_1d_char_array(int order, const char *name) {

       char *array = (char *) smalloc(order*sizeof(char),name);
       return array;

}
////////////////////////////////////////////////////////////////////////
/* Deallocation of a 1d character array */
void destroy_1d_char_array(char *array) {

     if (array == NULL) {
        printf("In destroy_1d_char_array pointer is NULL\n");
        return;
     }
     sfree(array);

}
////////////////////////////////////////////////////////////////////////
/* Dynamic allocation of a 2d character array */
char **create_2d_char_array(int n1, int n2, const char *name) {

       char *data = (char *) smalloc(n1*n2*sizeof(char),name);
       char **array = (char **) smalloc(n1*sizeof(char *),name);

       int order = 0, i;
       for (i = 0; i < n1; i++) {
           array[i] = &data[order];
           order += n2;
       }
       return array;

}
////////////////////////////////////////////////////////////////////////
/* Deallocation of a 2d character array */
void destroy_2d_char_array(char **array) {

     if (array == NULL) {
        printf("In destroy_2d_char_array pointer is NULL\n");
        return;
     }
     sfree(array[0]);
     sfree(array);

}

