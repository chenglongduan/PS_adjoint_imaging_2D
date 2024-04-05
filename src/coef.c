/* only provide 8th order coefs, details depend on Max error in group velocity */

#include "funclist.h"

float coef(float *hc, int maxerror, int fdorder, int MYID)
{
    /* Format: {min grids per shortest wavelength, 4 coefficients corresponding to L=4 (8th-order FD)} */
    float taylor_8th[7]      ={6.0, 1225.0/1024.0,     -245.0/3072.0,     49.0/5120.0,       -5.0/7168.0,       0.0,             0.0}; /*Taylor*/
    float holberg_8th_01[7]  ={3.68,1.225631165324630,-0.099498903588899,0.018041775664497,-0.002619544946100, 0.0, 0.0}; /*Holberg, 0.1% group velocity error*/
    float holberg_8th_002[7] ={4.41,1.217079526732709,-0.093278813464542,0.014936906894527,-0.001730650408702, 0.0, 0.0}; /*Holberg, 0.02% group velocity error*/

    
    float sum=0.0;
    int i;
    
    if(fdorder!=8 && MYID==0) {printf("fdorder can only be 8!");MPI_Finalize();exit(1);}

    switch(maxerror){
        case 0:
        {
            for (i=1;i<7;i++){
                sum += fabs(taylor_8th[i]);
            }
            memcpy(hc,taylor_8th,7*sizeof(float));
            break;
        }
        case 1:
        {
            for (i=1;i<7;i++){
                sum += fabs(holberg_8th_01[i]);
            }
            memcpy(hc,holberg_8th_01,7*sizeof(float));
            break;
        }
        case 2:
        {
            for (i=1;i<7;i++){
                sum += fabs(holberg_8th_002[i]);
            }
            memcpy(hc,holberg_8th_002,7*sizeof(float));
            break;
        }
    }
      
    
    
    return sum;
}
