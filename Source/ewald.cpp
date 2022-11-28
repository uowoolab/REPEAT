#include "ewald.h"
using namespace std;
////////////////////////////////////////////////////////////////////////////
/* Initializing Ewald data sets */
double init_ewald_data(void) {

     int i_dim, j_dim;
     double dist, alpha;

     if(flag_cutoff == 0) R_cutoff = 0.5*pow(Box_Volume,1/3);
     printf("Cutoff radius = %f \n",R_cutoff); 

     for(j_dim=0;j_dim<NDIM;j_dim++){    
        dist = 0.0;
        for(i_dim=0;i_dim<NDIM;i_dim++){
           dist += pow(real_box_vector[j_dim][i_dim],2);
        }
        dist = sqrt(dist);
        NMAX[j_dim] = static_cast<int> (floor(R_cutoff/dist)+1);
        //printf("# of cells in %d direction = %d\n",j_dim+1,NMAX[j_dim]);
     }
     alpha = pow(PI/R_cutoff,2);
     printf("Width of Gaussian in Ewald summation (alpha) = %f\n",alpha);
     k_vectors_coeff(alpha);
     return alpha;    
 
}
///////////////////////////////////////////////////////////////
/* Creating the array of k vectors */
void k_vectors_coeff(double alpha) {

     int TOTK, KX, KY, KZ, KSQ;
     double beta, RKX, RKY, RKZ, RKSQ;

     // Memory allocation
     KVEC = create_1d_double_array(MAXK,"k vectors");

     beta = 1.0/(4*alpha);
     TOTK = -1;
     for(KX=0;KX<=KMAX;KX++) {
        for(KY=-KMAX;KY<=KMAX;KY++) {
           for(KZ=-KMAX;KZ<=KMAX;KZ++) {
              RKX = recip_box_vector[0][0]*KX + \
              recip_box_vector[1][0]*KY + \
              recip_box_vector[2][0]*KZ;

              RKY = recip_box_vector[0][1]*KX + \
              recip_box_vector[1][1]*KY + \
              recip_box_vector[2][1]*KZ;

              RKZ = recip_box_vector[0][2]*KX + \
              recip_box_vector[1][2]*KY + \
              recip_box_vector[2][2]*KZ;

              KSQ = KX*KX + KY*KY + KZ*KZ;
              if((KSQ < KSQMAX) && (KSQ != 0)) {
                TOTK++;
                if(TOTK > MAXK) {
                  printf("ERROR: KVEC IS TOO SMALL\n");
                  exit(EXIT_FAILURE);
                }
                RKSQ = RKX*RKX + RKY*RKY + RKZ*RKZ;
                KVEC[TOTK] = (4*PI/Box_Volume)*exp(-beta*RKSQ)/RKSQ;
                //printf(" %le\n",KVEC[TOTK]);
              }
           }
        }
     }

}
///////////////////////////////////////////////////////////////////////
/* Computing real term of the Ewald sum */
double Phi_real_coeff(double alpha,double *pos_grid,double *pos_charge) {

       int i_cell_x, i_cell_y, i_cell_z, i_dim;
       double dist, Phi, delta_dist[NDIM];

       Phi = 0.0;
#pragma ivdep
       for(i_cell_x=-NMAX[0];i_cell_x<=NMAX[0];i_cell_x++) {
#pragma ivdep
          for(i_cell_y=-NMAX[1];i_cell_y<=NMAX[1];i_cell_y++) {
#pragma ivdep
             for(i_cell_z=-NMAX[2];i_cell_z<=NMAX[2];i_cell_z++) {
                dist = 0.0;
                for(i_dim=0;i_dim<NDIM;i_dim++) {
                   delta_dist[i_dim] = pos_grid[i_dim] - \
                   (pos_charge[i_dim] + \
                   i_cell_x*real_box_vector[0][i_dim] + \
                   i_cell_y*real_box_vector[1][i_dim] + \
                   i_cell_z*real_box_vector[2][i_dim]);
                   dist += pow(delta_dist[i_dim],2);
                }
                dist = sqrt(dist);
                if(dist <= R_cutoff) {
                   Phi += erfc(sqrt(alpha)*dist)/dist;
                }
             }
          }
       }
       return Phi;

}
////////////////////////////////////////////////////////////////////////
/* Computing the reciprocal space term of the Ewald sum */
double NEW_Phi_recp_coeff(double *delta_r) {

       int TOTK, KX, KY, KZ, KSQ, index_1, index_2;
       double RX, RY, RZ, FACTOR, Phi;
       complex<double> EIKX[KMAX+1], EIKY[2*KMAX+1], EIKZ[2*KMAX+1], \
       EIKR, VD;
       
       RX = delta_r[0];
       RY = delta_r[1];
       RZ = delta_r[2];      

       EIKX[0] = complex<double>(1.0,0.0);
       index_1 = relabel_index(0);
       EIKY[index_1] = complex<double>(1.0,0.0);
       EIKZ[index_1] = complex<double>(1.0,0.0);

       EIKX[1] = complex<double>(cos(recip_box_vector[0][0]*RX + \
       recip_box_vector[0][1]*RY + recip_box_vector[0][2]*RZ), \
       sin(recip_box_vector[0][0]*RX + \
       recip_box_vector[0][1]*RY + recip_box_vector[0][2]*RZ));

       index_1 = relabel_index(1);
       EIKY[index_1] = complex<double>(cos(recip_box_vector[1][0]*RX + \
       recip_box_vector[1][1]*RY + recip_box_vector[1][2]*RZ), \
       sin(recip_box_vector[1][0]*RX + \
       recip_box_vector[1][1]*RY + recip_box_vector[1][2]*RZ));

       EIKZ[index_1] = complex<double>(cos(recip_box_vector[2][0]*RX + \
       recip_box_vector[2][1]*RY + recip_box_vector[2][2]*RZ), \
       sin(recip_box_vector[2][0]*RX + \
       recip_box_vector[2][1]*RY + recip_box_vector[2][2]*RZ));

       index_2 = relabel_index(-1);
       EIKY[index_2] = conj(EIKY[index_1]);
       EIKZ[index_2] = conj(EIKZ[index_1]);

       for(KX=2;KX<=KMAX;KX++) {
          EIKX[KX] =  EIKX[KX-1]*EIKX[1];
       } 
       
       for(KY=2;KY<=KMAX;KY++) {
          index_1 = relabel_index(KY);
          EIKY[index_1] = EIKY[relabel_index(KY-1)]*EIKY[relabel_index(1)];
          EIKY[relabel_index(-KY)] = conj(EIKY[index_1]);
       }

       for(KZ=2;KZ<=KMAX;KZ++) {
          index_1 = relabel_index(KZ);
          EIKZ[index_1] = EIKZ[relabel_index(KZ-1)]*EIKZ[relabel_index(1)];
          EIKZ[relabel_index(-KZ)] = conj(EIKZ[index_1]);
       }

       VD = complex<double>(0.0,0.0);
       TOTK = -1;
       FACTOR = 1.0;
       KX=0;
          for(KY=-KMAX;KY<=KMAX;KY++) {
             index_1 = relabel_index(KY);
             for(KZ=-KMAX;KZ<=KMAX;KZ++) {
                index_2 = relabel_index(KZ);
                KSQ = KX*KX + KY*KY + KZ*KZ;
                if((KSQ < KSQMAX) && (KSQ != 0)) {
                  TOTK++;
                  EIKR = EIKX[KX]*EIKY[index_1]*EIKZ[index_2];     
                  VD += FACTOR*KVEC[TOTK]*EIKR;
                }
             }
          }
       FACTOR = 2.0;
#pragma ivdep
       for(KY=-KMAX;KY<=KMAX;KY++) {
          V_my_index_1[KY] = relabel_index(KY);
       }
#pragma ivdep
       for(KZ=-KMAX;KZ<=KMAX;KZ++) {
                V_my_index_2[KZ] = relabel_index(KZ);
       }
#pragma ivdep
       for(KX=1;KX<=KMAX;KX++) {
#pragma ivdep
          for(KY=-KMAX;KY<=KMAX;KY++) {
#pragma ivdep
             for(KZ=-KMAX;KZ<=KMAX;KZ++) {
                KSQ = KX*KX + KY*KY + KZ*KZ;
                if((KSQ < KSQMAX) && (KSQ != 0)) {
                  TOTK++;
                  EIKR = EIKX[KX]*EIKY[V_my_index_1[KY]]*EIKZ[V_my_index_2[KZ]];     
                  VD += FACTOR*KVEC[TOTK]*EIKR;
                }
             }
          }
       }
       Phi = VD.real();
       return Phi;       

}
double Phi_recp_coeff(double *delta_r) {

       int TOTK, KX, KY, KZ, KSQ, index_1, index_2;
       double RX, RY, RZ, FACTOR, Phi;
       complex<double> EIKX[KMAX+1], EIKY[2*KMAX+1], EIKZ[2*KMAX+1], \
       EIKR, VD;
       
       RX = delta_r[0];
       RY = delta_r[1];
       RZ = delta_r[2];      

       EIKX[0] = complex<double>(1.0,0.0);
       index_1 = relabel_index(0);
       EIKY[index_1] = complex<double>(1.0,0.0);
       EIKZ[index_1] = complex<double>(1.0,0.0);

       EIKX[1] = complex<double>(cos(recip_box_vector[0][0]*RX + \
       recip_box_vector[0][1]*RY + recip_box_vector[0][2]*RZ), \
       sin(recip_box_vector[0][0]*RX + \
       recip_box_vector[0][1]*RY + recip_box_vector[0][2]*RZ));

       index_1 = relabel_index(1);
       EIKY[index_1] = complex<double>(cos(recip_box_vector[1][0]*RX + \
       recip_box_vector[1][1]*RY + recip_box_vector[1][2]*RZ), \
       sin(recip_box_vector[1][0]*RX + \
       recip_box_vector[1][1]*RY + recip_box_vector[1][2]*RZ));

       EIKZ[index_1] = complex<double>(cos(recip_box_vector[2][0]*RX + \
       recip_box_vector[2][1]*RY + recip_box_vector[2][2]*RZ), \
       sin(recip_box_vector[2][0]*RX + \
       recip_box_vector[2][1]*RY + recip_box_vector[2][2]*RZ));

       index_2 = relabel_index(-1);
       EIKY[index_2] = conj(EIKY[index_1]);
       EIKZ[index_2] = conj(EIKZ[index_1]);

       #pragma ivdep
       for(KX=2;KX<=KMAX;KX++) {
          EIKX[KX] =  EIKX[KX-1]*EIKX[1];
       } 
       
       #pragma ivdep
       for(KY=2;KY<=KMAX;KY++) {
          index_1 = relabel_index(KY);
          EIKY[index_1] = EIKY[relabel_index(KY-1)]*EIKY[relabel_index(1)];
          EIKY[relabel_index(-KY)] = conj(EIKY[index_1]);
       }

       #pragma ivdep
       for(KZ=2;KZ<=KMAX;KZ++) {
          index_1 = relabel_index(KZ);
          EIKZ[index_1] = EIKZ[relabel_index(KZ-1)]*EIKZ[relabel_index(1)];
          EIKZ[relabel_index(-KZ)] = conj(EIKZ[index_1]);
       }

       VD = complex<double>(0.0,0.0);
       TOTK = -1;
       #pragma ivdep
       for(KX=0;KX<=KMAX;KX++) {
          if(KX == 0) {
            FACTOR = 1.0;
          } else {
            FACTOR = 2.0;
          }
          #pragma ivdep
          for(KY=-KMAX;KY<=KMAX;KY++) {
             index_1 = relabel_index(KY);
          #pragma ivdep
             for(KZ=-KMAX;KZ<=KMAX;KZ++) {
                index_2 = relabel_index(KZ);
                KSQ = KX*KX + KY*KY + KZ*KZ;
                if((KSQ < KSQMAX) && (KSQ != 0)) {
                  TOTK++;
                  EIKR = EIKX[KX]*EIKY[index_1]*EIKZ[index_2];     
                  VD += FACTOR*KVEC[TOTK]*EIKR;
                }
             }
          }
       }
       Phi = VD.real();
       return Phi;       

}
////////////////////////////////////////////////////////////////////////
/* Computing Coulomb term for molecular systems */
double Phi_coulomb(double *pos_grid,double *pos_charge) {

       int i_dim;
       double Phi, dist;

       dist = 0.0;
       for(i_dim=0;i_dim<NDIM;i_dim++) {
          dist += pow((pos_grid[i_dim] - pos_charge[i_dim]),2);
       }
       dist = sqrt(dist);
       Phi = 1.0/dist;
       return Phi;

}

