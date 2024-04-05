


#include "funclist.h"


void pad_bound ( int nx, int nz, int nxpad, int nzpad, int nb, float **raw_mod, float **pad_mod )

{
    
    int ix, iz;
    
    for     (iz=0; iz<nz; iz++){
        for (ix=0; ix<nx; ix++){
            pad_mod[iz+nb][ix+nb] = raw_mod[iz][ix];  /* on normal region */
        }
    }
    
    for     (iz=0; iz<nzpad; iz++){
        for (ix=0; ix<nb;    ix++){
            pad_mod[iz][      ix  ] = pad_mod[iz][      nb  ];    /* pad left */
            pad_mod[iz][nxpad-ix-1] = pad_mod[iz][nxpad-nb-1];    /* pad right */
        }
    }
    
    for     (iz=0; iz<nb;    iz++){
        for (ix=0; ix<nxpad; ix++){
            pad_mod[      iz  ][ix] = pad_mod[      nb  ][ix];      /* pad top */
            pad_mod[nzpad-iz-1][ix] = pad_mod[nzpad-nb-1][ix];      /* pad bottom */
        }
    }



}
