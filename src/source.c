#include "funclist.h"


int **source(int *nsrc, char *SOURCE_FILE, int MYID, float dx, float dz, int nx, int nz, int nb)
{

    FILE *fpsrc;
    int isrc, k;
    char cline[256];
    float xsrc, zsrc;
    int ixsrc, izsrc;
    int **srcpos=NULL;
    
    
    if (MYID==0){
        printf("\n Reading source file: %s\n", SOURCE_FILE);
        
        fpsrc = fopen(SOURCE_FILE, "r");
        if (fpsrc==NULL) {printf(" Source file could not be opened !");MPI_Finalize();exit(1);}
        
        *nsrc = 0;
        fscanf(fpsrc,"%d",nsrc);
        printf(" Number of sources: %d \n",*nsrc);
        
        srcpos = sf_intalloc2(*nsrc,2); memset(srcpos[0],0,2*(*nsrc)*sizeof(int));
        
        rewind(fpsrc);
        fgets(cline, 255, fpsrc);
        
        for (isrc=0; isrc<(*nsrc); isrc++){
            fgets(cline, 255, fpsrc);
            sscanf(cline, "%f%f", &xsrc, &zsrc); //real physical location (do not consider boundary and C convention)
            ixsrc = round(xsrc/dx) + 1;
            izsrc = round(zsrc/dz) + 1;
            if (ixsrc>nx || izsrc>nz){
                sf_warning("Err: Sources exceed the computing zone!");
            }
            srcpos[0][isrc] = ixsrc + nb - 1; //ix
            srcpos[1][isrc] = izsrc + nb - 1; //iz
        }
        fclose(fpsrc); //---
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_Bcast(nsrc, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    if (MYID!=0) {srcpos = sf_intalloc2(*nsrc,2); memset(srcpos[0],0,2*(*nsrc)*sizeof(int));}
    MPI_Bcast(&srcpos[0][0], (*nsrc)*2, MPI_INT, 0, MPI_COMM_WORLD);
    
    if(MYID==0){
        sf_warning(" Source list (C conv w/ nb):\n");
        printf(" ix  \tiz \n");
        printf(" --  \t-- \n");
        if(*nsrc>25){
            for (k=0;k<12;k++)
                printf(" %d \t %d\n",srcpos[0][k],srcpos[1][k]);
            printf(" ... \t ... \n");
            for (k=*nsrc-12;k<*nsrc;k++)
                printf(" %d \t %d\n",srcpos[0][k],srcpos[1][k]);
            printf("\n\n");
	}else{
	    for (k=0;k<*nsrc;k++){
                printf(" %d \t %d\n",srcpos[0][k],srcpos[1][k]);
	    }
	    printf("\n\n");
	}
    }
    
    
    return srcpos;

}
