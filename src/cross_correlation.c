
#include "funclist.h"



void cross_correlation( float **num, float **den, float **sp, float **gp, int nb, int fdo, int nx, int nz )
/*< cross correlation >*/
{
    
    int ix, iz;
    float xcorr;

    for(iz=0; iz<nz; iz++){
	for(ix=0; ix<nx; ix++){
            
            /*cross-correlation at it-step*/
	    xcorr = sp[iz+nb+fdo][ix+nb+fdo] * gp[iz+nb+fdo][ix+nb+fdo];
	    
	    /*time integral == incident angle stacking*/
	    num[iz][ix] += xcorr;
	    //den[iz][ix] += gp[iz+nb+fdo][ix+nb+fdo] * gp[iz+nb+fdo][ix+nb+fdo]; //receiver-side illumination
            den[iz][ix] += sp[iz+nb+fdo][ix+nb+fdo] * sp[iz+nb+fdo][ix+nb+fdo]; //source-side illumination
            //den[iz][ix] += sp[iz+nb+fdo][ix+nb+fdo] * sp[iz+nb+fdo][ix+nb+fdo] * gp[iz+nb+fdo][ix+nb+fdo] * gp[iz+nb+fdo][ix+nb+fdo]; //diagonal Hessian
	    
	}
    }

}
