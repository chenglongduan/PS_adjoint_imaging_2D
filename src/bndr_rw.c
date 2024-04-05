#include "funclist.h"


void bndr_rw(bool read, float **p, char *tmp_srcwfd_dir, int it, int isrc, int fdo, int nb, int nx, int nz)
{
    
    int iz, ix;
    FILE *fr=NULL, *fw=NULL;
    char filename[512];
    float value;

    if(read){ //read in

	sprintf(filename,"%sT%d.S%d",tmp_srcwfd_dir,it,isrc); //follow the loop count, so not +1
	fr=fopen(filename,"r");
	if (fr==NULL) {sf_warning(" Could not open source wavefield file! ");MPI_Finalize();exit(1);}
	
	for (ix=fdo+nb; ix<nx+nb+fdo; ix++){
	    for (iz=fdo+nb; iz<nz+nb+fdo; iz++){
	        fread(&value, sizeof(float), 1, fr);				
		if (isnan(value)) {
           	    sf_warning(" Found NaN-Values in the saved source wavefield!");
           	    MPI_Finalize();exit(1);
           	}			
		p[iz][ix]=value;
	    }
	}
	fclose(fr);
	
    }
    else{  //write out	
	
	sprintf(filename,"%sT%d.S%d",tmp_srcwfd_dir,it,isrc); //follow the loop count, so not +1
        fw = fopen(filename,"wb");
        if (fw==NULL) {sf_warning(" Could not write source wavefield file! ");MPI_Finalize();exit(1);}
        
        for (ix=fdo+nb; ix<nx+nb+fdo; ix++){
            for (iz=fdo+nb; iz<nz+nb+fdo; iz++){
                fwrite(&p[iz][ix],sizeof(float),1,fw);
            }
        }
        fclose(fw);

    }

}


