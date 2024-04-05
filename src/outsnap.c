
#include "funclist.h"


void outsnap( int snaptype, float **p, float **vx, float **vz, int isrc, int it, int nxpad, int nzpad, int nb, int fdo, char *snapdir, char ext[50] )

{
    
    FILE *fp_snap=NULL;
    char name_buffer[512];
    int ix, iz;
    
    
    switch (snaptype) {
        case 1: //pressure wavefield
            sprintf(name_buffer,"%s%s%d%s%d%s",snapdir,"p",isrc+1,"snap",it+1,ext);
            fp_snap = fopen(name_buffer,"wb");
            for (ix=fdo+nb; ix<nxpad+fdo-nb; ix++){
                for (iz=fdo+nb; iz<nzpad+fdo-nb; iz++){
                    fwrite(&p[iz][ix],sizeof(float),1,fp_snap);
                }
            }
            fclose(fp_snap);
            break;
        case 2: //vx wavefield
            sprintf(name_buffer,"%s%s%d%s%d%s",snapdir,"vx",isrc+1,"snap",it+1,ext);
            fp_snap = fopen(name_buffer,"wb");
            for (ix=fdo+nb; ix<nxpad+fdo-nb; ix++){
                for (iz=fdo+nb; iz<nzpad+fdo-nb; iz++){
                    fwrite(&vx[iz][ix],sizeof(float),1,fp_snap);
                }
            }
            fclose(fp_snap);
            break;
        case 3: //vz wavefield
            sprintf(name_buffer,"%s%s%d%s%d%s",snapdir,"vz",isrc+1,"snap",it+1,ext);
            fp_snap = fopen(name_buffer,"wb");
            for (ix=fdo+nb; ix<nxpad+fdo-nb; ix++){
                for (iz=fdo+nb; iz<nzpad+fdo-nb; iz++){
                    fwrite(&vz[iz][ix],sizeof(float),1,fp_snap);
                }
            }
            fclose(fp_snap);
            break;
        
    }
    
    
    
    
    


}
