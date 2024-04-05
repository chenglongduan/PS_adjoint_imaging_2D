
#include "funclist.h"


void outsnapimg(float **num, float **den, int it, int MYID, int nx, int nz, char *snapimg_dir) //clduan
{

    FILE *fp_snap=NULL;
    char name_buffer[512];
    int ix, iz;
    float **img_1src_it=NULL, **img_nsrc_it=NULL; //clduan
    
    
    sprintf(name_buffer,"%simg.t%d",snapimg_dir,it+1); //clduan
    
    /*allocation*/
    img_1src_it=sf_floatalloc2(nx,nz); memset(img_1src_it[0], 0, nz*nx*sizeof(float));
    img_nsrc_it=sf_floatalloc2(nx,nz); memset(img_nsrc_it[0], 0, nz*nx*sizeof(float));
    
    /*image at it-time & is-shot*/
    for(iz=0; iz<nz; iz++){
        for(ix=0; ix<nx; ix++){
            img_1src_it[iz][ix] = num[iz][ix]/(den[iz][ix]+SF_EPS);
        }
    }
    
    /*local stacking*/
    MPI_Allreduce(&img_1src_it[0][0], &img_nsrc_it[0][0], nx*nz, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    
    /*output*/
    if(MYID==0){
        fp_snap = fopen(name_buffer,"wb");
        for(ix=0; ix<nx; ix++){
            for(iz=0; iz<nz; iz++){
                fwrite(&img_nsrc_it[iz][ix],sizeof(float),1,fp_snap);
            }
        }
        fclose(fp_snap);
    }
    
    /*free matrix*/
    free(*img_1src_it);free(img_1src_it);
    free(*img_nsrc_it);free(img_nsrc_it);
    
    /*The following is the previous version*/
    /*sprintf(name_buffer,"%simg.s%d.t%d",snapimg_dir,isrc+1,it+1);
    fp_snap = fopen(name_buffer,"wb");
    
    for(ix=0; ix<nx; ix++){
        for(iz=0; iz<nz; iz++){
            
            img_1src_it[iz][ix] = num[iz][ix]/(den[iz][ix]+SF_EPS);
            
            fwrite(&img_1src_it[iz][ix],sizeof(float),1,fp_snap);
        
        }
    }
    fclose(fp_snap);*/
    

}
