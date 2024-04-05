#include "funclist.h"


int **receiver(int *ntr, char *REC_FILE, int MYID, float dx, float dz, int nx, int nz, int nb)
{

    FILE *fpr;
    int itr=1, k;
    char bufferstring[10], buffer[512];
    bool testbuff1, testbuff2, testbuff3;
    float xrec, zrec;
    int ixrec, izrec;
    int **recpos=NULL;
    
    
    if(MYID==0)  printf("\n Process %d is reading receiver file %s\n",MYID,REC_FILE);
        
    fpr = fopen(REC_FILE, "r");
    if (fpr==NULL) {printf("Receiver file could not be opened!");MPI_Finalize();exit(1);}
        
    /* count the number of receivers */
    *ntr = 0;
    while(fgets(buffer, 512, fpr)){
        testbuff1 = strchr(buffer,'#');
        testbuff2 = strchr(buffer,'%');
        testbuff3 = sscanf(buffer,"%s",bufferstring)==1;
        /* checks if the line contains '%' or '#' which is a comment line; and if reading a string is successful, which is not an empty line */
        if (((testbuff1==1 || testbuff2==1)==0) && testbuff3==1) ++(*ntr);
    }
    rewind(fpr);
        
    recpos = sf_intalloc2(*ntr,2); memset(recpos[0],0,2*(*ntr)*sizeof(int));
        
    for (itr=0; itr<(*ntr); itr++){
        fscanf(fpr, "%f%f\n", &xrec, &zrec); //real physical location (do not consider boundary and C convention)
            
        ixrec = round(xrec/dx) + 1;
        izrec = round(zrec/dz) + 1;
        if (ixrec>nx || izrec>nz){
            sf_warning("Err: Receivers exceed the computing zone!");
        }
            
        recpos[0][itr] = ixrec + nb - 1; //ix
        recpos[1][itr] = izrec + nb - 1; //iz
    }
    fclose(fpr);
        
    if(MYID==0)  printf(" Process %d: number of receivers = %i\n",MYID,*ntr);
    
        
    if(MYID==0){
        sf_warning(" Process %d: receiver list (C conv w/ nb):\n",MYID);
        printf(" ix  \tiz \n");
        printf(" --  \t-- \n");
        if(*ntr>25){
            for (k=0;k<12;k++)
                printf(" %d \t %d\n",recpos[0][k],recpos[1][k]);
            printf(" ... \t ... \t \n");
            for (k=*ntr-12;k<*ntr;k++)
                printf(" %d \t %d\n",recpos[0][k],recpos[1][k]);
            printf("\n");
        }else{
            for (k=0;k<*ntr;k++){
                printf(" %d \t %d\n",recpos[0][k],recpos[1][k]);
            }
            printf("\n");
        }

    }

    return recpos;

}
