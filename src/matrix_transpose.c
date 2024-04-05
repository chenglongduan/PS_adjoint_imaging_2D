
#include "funclist.h"


//For input (in) and output (ex) translation.
//Both input and output are based on the SAME coordinate system. Thus, nx,ny,nz are unique

/*coordinate:
       _ _ _ _ X
      |
      |
      |
      Z           
*/
void matrix_transp_2D(float **mat_ex, float **mat_in, int nx, int nz, int io)
{
    int ix, iz;
    if(io==1){ //ex->in
        for(ix=0; ix<nx; ix++){
	    for(iz=0; iz<nz; iz++)
	    {
	        mat_in[iz][ix]=mat_ex[ix][iz];
	    }
        }
    }
    if(io==2){ //in->ex
        for(iz=0; iz<nz; iz++){
	    for(ix=0; ix<nx; ix++){
                mat_ex[ix][iz]=mat_in[iz][ix];
            }
        }
    }
}


