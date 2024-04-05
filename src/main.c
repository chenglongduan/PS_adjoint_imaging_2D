/* 2D Acoustic RTM (adjoint imaging) with CPML boundary */
/* First-order stagger-grid FD, 8th space + 2nd time */

#include "funclist.h"
    


int main(int argc, char* argv[])
{
    
    int nb, fdorder, pstp, maxerror,ng,snaptime,snapimg,snaptype,seismo_cal;
    int nz,nx, fdo, nxpad,nzpad,nt,isrc,nsrc,isx,isz,igx,igz;
    int i,j,src_func,outsrc,it,irec;
    int iz,ix,c1,c3,c4,do_srcwfd,out_img_byshot;
    int **srcpos=NULL, **recpos=NULL;
    float *wlt=NULL, *hc=NULL;
    float **vv1=NULL,**vv2=NULL,**rr=NULL,**velfwd=NULL,**veladj=NULL,**rho=NULL,**velfwd_pad=NULL,
    **veladj_pad=NULL,**rhopad=NULL,**inv_rho_half_x=NULL,**inv_rho_half_z=NULL;
    char *wlt_name=NULL,*outdir_src=NULL,*snapdir=NULL,*snapimg_dir=NULL,*dat_dir_prefix=NULL,*cal_dir_prefix=NULL,*imgName=NULL,*tmp_srcwfd_dir=NULL,
    *source_file=NULL,*recv_dir=NULL,*img_byshot_dir=NULL;
    float f0, VPPML,oz,ox,dz,dx,amp,dt,tt,td=0.0,tau=0.0,gamma,vmax,vmin,c2,c5;
    sf_file fvel1, fvel2, fden, fwlt, fdat, frtm;
    char name_buffer[512],datName[512],calName[512],receiver_file[512],tmp_img1[512];
    FILE *fp,*fdcal,*fimg1;
    float *d_x=NULL, *K_x=NULL, *alpha_prime_x=NULL, *a_x=NULL, *b_x=NULL, *d_z=NULL, *K_z=NULL, *alpha_prime_z=NULL, *a_z=NULL, *b_z=NULL, 
    *d_x_half=NULL, *K_x_half=NULL, *alpha_prime_x_half=NULL, *a_x_half=NULL, *b_x_half=NULL, 
    *d_z_half=NULL, *K_z_half=NULL, *alpha_prime_z_half=NULL, *a_z_half=NULL, *b_z_half=NULL;
    float **sp=NULL, **svx=NULL, **svz=NULL, **gp=NULL, **gvx=NULL, **gvz=NULL;
    float **psi1_dp_dx=NULL, **psi1_dp_dz=NULL, **psi1_dvx_dx=NULL, **psi1_dvz_dz=NULL, **psi2_dp_dx=NULL, **psi2_dp_dz=NULL, **psi2_dvx_dx=NULL, **psi2_dvz_dz=NULL;
    float **dobs=NULL, **dcal=NULL, **den=NULL, **num=NULL, **img_1src=NULL, **img_nsrc=NULL, **img_out=NULL;

    /* MPI variables */
    int NP,MYID;
    
    
    /* Initialize MPI environment: 
    distribute main() to NP cores; each core is assigned an ID from 0 to NP-1 */
    MPI_Init ( &argc, &argv );
    MPI_Comm_size ( MPI_COMM_WORLD, &NP );
    MPI_Comm_rank ( MPI_COMM_WORLD, &MYID ); //0:NP-1
    
    /* initialize Madagascar */
    sf_init(argc,argv);


    /* CPML parameters */
    if(!sf_getint("nb", &nb)) nb = 30;
    if(!sf_getfloat("f0", &f0)) f0 = 30.0f;
    if(!sf_getfloat("VPPML", &VPPML)) VPPML = 4100.0f;
    
    /* Finite difference */
    if(!sf_getint("fdorder", &fdorder)) fdorder = 8;
    if(!sf_getint("printstep", &pstp)) pstp = 400;
    if(!sf_getint("maxerror", &maxerror)) maxerror = 0;
    
    /* Model space */
    fvel1 = sf_input("velocity_fwd"); /*un-padded forward model*/
    fvel2 = sf_input("velocity_adj"); /*un-padded adjoint model*/
    fden = sf_input("density");  /*un-padded model*/
    if(!sf_histint  (fvel1, "n1", &nz)) {fprintf(stderr,"No n1= in velocity\n");MPI_Finalize();exit(1);}
    if(!sf_histint  (fvel1, "n2", &nx)) {fprintf(stderr,"No n2= in velocity\n");MPI_Finalize();exit(1);}
    if(!sf_histfloat(fvel1, "o1", &oz)) {fprintf(stderr,"No o1= in velocity\n");MPI_Finalize();exit(1);}
    if(!sf_histfloat(fvel1, "o2", &ox)) {fprintf(stderr,"No o2= in velocity\n");MPI_Finalize();exit(1);}
    if(!sf_histfloat(fvel1, "d1", &dz)) {fprintf(stderr,"No d1= in velocity\n");MPI_Finalize();exit(1);}
    if(!sf_histfloat(fvel1, "d2", &dx)) {fprintf(stderr,"No d2= in velocity\n");MPI_Finalize();exit(1);}
    
    fdo = fdorder/2;
    nxpad = nx+2*nb;
    nzpad = nz+2*nb;
    
    /* Time step */
    if(!sf_getfloat("dt", &dt)) dt = 3.5e-4;
    if(!sf_getfloat("tt", &tt)) tt = 1.0f;  /*total modeling time*/
    
    nt = (int)ceil(tt/dt);
    
    /* =========================================Source=========================================== */
    /* source locations */
    source_file = sf_charalloc(512);
    source_file = sf_getstring("source_file");
    srcpos = source(&nsrc, source_file, MYID, dx, dz, nx, nz, nb);
    
    /* source wavelet */
    if(!sf_getfloat("amp", &amp)) amp = 1.e4; /* amplitude scalar of source wavelet*/
    if(!sf_getint("src_func", &src_func)) src_func = 1; /*Assign source time function*/
    wlt = sf_floatalloc(nt);
    memset(&wlt[0],0,nt*sizeof(float)); 
    if(src_func==3){
        if(MYID==0) sf_warning("STF: Customize\n");
        wlt_name = sf_charalloc(512);
        if( NULL== (wlt_name=sf_getstring("src_wlt")) ){ 
          fprintf(stderr,"No external source wavelet inputs"); 
          MPI_Finalize(); exit(1); 
        }
        fwlt = sf_input(wlt_name);
        if(!sf_histint  (fwlt, "n1", &c1)) {fprintf(stderr,"No n1= in wavelet\n");MPI_Finalize();exit(1);}
        if(!sf_histfloat(fwlt, "d1", &c2)) {fprintf(stderr,"No d1= in wavelet\n");MPI_Finalize();exit(1);}
        if(c1!=nt || c2!=dt) {printf("wavelet nt,dt wrong");MPI_Finalize();exit(1);}
        sf_floatread(wlt,nt,fwlt);
    }else{
        if(src_func==1 && MYID==0) sf_warning("STF: Ricker\n");
        if(src_func==2 && MYID==0) sf_warning("STF: Gaussian\n");
        if(!sf_getfloat("time_delay", &td)) td = 0.0;
    }
    for (it=0; it<nt; it++){            
        switch (src_func){
              case 1: /* Ricker */
                tau = PI*((float)it*dt-1.5/f0-td)*f0;
                wlt[it] = ((1.0-2.0*tau*tau)*exp(-tau*tau)) * amp;
              break;
              case 2: /* Gaussian */
                tau = PI*((float)it*dt-1.5/f0-td)*f0;
                wlt[it] = exp(-tau*tau)/(2.0*PI*PI*f0*f0) * amp;
              break;
              case 3: /* Customize */
                wlt[it] = wlt[it] * amp;
              break;
        }
    }
    /*store source wavelet*/
    if(!sf_getint("outsrc", &outsrc)) outsrc = 0;
    if(outsrc==1 && MYID==0){
        outdir_src = sf_charalloc(512);
        outdir_src = sf_getstring("outdir_src");
        sprintf(name_buffer,"%s%s%s",outdir_src,"src_wlt",".bin");
        fp = fopen(name_buffer,"wb");
	for(it=0;it<nt;it++){
            fwrite(&wlt[it],sizeof(float),1,fp);
        }
       
        fclose(fp);
    }
    /* time integration for true amplitude RTM */
    //for(it=1;it<nt;it++) wlt[it]=wlt[it]+wlt[it-1]; //remove by clduan
    
    
    /* =========================================Receiver=========================================== */
    recv_dir = sf_charalloc(512);
    recv_dir = sf_getstring("recv_dir");
    
    /* =========================================Models=========================================== */
    vv1 = sf_floatalloc2(nz, nx); memset(vv1[0],0,nz*nx*sizeof(float));
    vv2 = sf_floatalloc2(nz, nx); memset(vv2[0],0,nz*nx*sizeof(float));
    rr = sf_floatalloc2(nz, nx); memset(rr[0],0,nz*nx*sizeof(float));
    
    velfwd = sf_floatalloc2(nx, nz); memset(velfwd[0],0,nz*nx*sizeof(float));
    veladj = sf_floatalloc2(nx, nz); memset(veladj[0],0,nz*nx*sizeof(float));
    rho = sf_floatalloc2(nx, nz); memset(rho[0],0,nz*nx*sizeof(float));
    
    velfwd_pad = sf_floatalloc2(nxpad, nzpad); memset(velfwd_pad[0],0,nzpad*nxpad*sizeof(float));
    veladj_pad = sf_floatalloc2(nxpad, nzpad); memset(veladj_pad[0],0,nzpad*nxpad*sizeof(float));
    rhopad = sf_floatalloc2(nxpad+1, nzpad+1); memset(rhopad[0],0,(nzpad+1)*(nxpad+1)*sizeof(float));   
    inv_rho_half_x = sf_floatalloc2(nxpad,nzpad); memset(inv_rho_half_x[0],0,nxpad*nzpad*sizeof(float));    
    inv_rho_half_z = sf_floatalloc2(nxpad,nzpad); memset(inv_rho_half_z[0],0,nxpad*nzpad*sizeof(float));
    
    // read
    sf_floatread(vv1[0],nz*nx,fvel1); /*forward velocity*/
    sf_floatread(vv2[0],nz*nx,fvel2); /*adjoint velocity*/
    sf_floatread(rr[0],nz*nx,fden); /*rho*/
    // transpose
    matrix_transp_2D(vv1, velfwd, nx, nz, 1);
    matrix_transp_2D(vv2, veladj, nx, nz, 1);
    matrix_transp_2D(rr, rho, nx, nz, 1);
    // pad velocity
    pad_bound ( nx, nz, nxpad, nzpad, nb, velfwd, velfwd_pad );
    pad_bound ( nx, nz, nxpad, nzpad, nb, veladj, veladj_pad );
    // pad density
    pad_bound ( nx, nz, nxpad, nzpad, nb, rho, rhopad );
    for (iz=0;iz<nzpad+1;iz++){
            rhopad[iz][nxpad] = rhopad[iz][nxpad-1]; /*1 more right boundary for interpolation*/
    }
    for (ix=0;ix<nxpad+1;ix++){
            rhopad[nzpad][ix] = rhopad[nzpad-1][ix]; /*1 more btm boundary for interpolation*/
    }    
    // 1/Rho_half
    for (iz=0;iz<nzpad;iz++){
        for (ix=0;ix<nxpad;ix++){	   
            inv_rho_half_z[iz][ix] = 1.0/(0.5 * (rhopad[iz][ix] + rhopad[iz+1][ix]));
            inv_rho_half_x[iz][ix] = 1.0/(0.5 * (rhopad[iz][ix] + rhopad[iz][ix+1]));
	}
    }
    
    /* =========================================CPML=========================================== */
    d_x = sf_floatalloc(nxpad); memset(d_x,0,nxpad*sizeof(float));    
    d_x_half = sf_floatalloc(nxpad); memset(d_x_half,0,nxpad*sizeof(float));    
    alpha_prime_x = sf_floatalloc(nxpad); memset(alpha_prime_x,0,nxpad*sizeof(float));    
    alpha_prime_x_half = sf_floatalloc(nxpad); memset(alpha_prime_x_half,0,nxpad*sizeof(float));   
    a_x = sf_floatalloc(nxpad); memset(a_x,0,nxpad*sizeof(float));   
    a_x_half = sf_floatalloc(nxpad); memset(a_x_half,0,nxpad*sizeof(float));   
    b_x = sf_floatalloc(nxpad); memset(b_x,0,nxpad*sizeof(float));   
    b_x_half = sf_floatalloc(nxpad); memset(b_x_half,0,nxpad*sizeof(float));
    K_x = sf_floatalloc(nxpad);
    K_x_half = sf_floatalloc(nxpad);
    for (i=0;i<nxpad;i++){
        K_x[i] = 1.0f;
        K_x_half[i] = 1.0f;
    }  
    d_z = sf_floatalloc(nzpad); memset(d_z,0,nzpad*sizeof(float));  
    d_z_half = sf_floatalloc(nzpad); memset(d_z_half,0,nzpad*sizeof(float));    
    alpha_prime_z = sf_floatalloc(nzpad); memset(alpha_prime_z,0,nzpad*sizeof(float));    
    alpha_prime_z_half = sf_floatalloc(nzpad); memset(alpha_prime_z_half,0,nzpad*sizeof(float));    
    a_z = sf_floatalloc(nzpad); memset(a_z,0,nzpad*sizeof(float));   
    a_z_half = sf_floatalloc(nzpad); memset(a_z_half,0,nzpad*sizeof(float));   
    b_z = sf_floatalloc(nzpad); memset(b_z,0,nzpad*sizeof(float));    
    b_z_half = sf_floatalloc(nzpad); memset(b_z_half,0,nzpad*sizeof(float));    
    K_z = sf_floatalloc(nzpad);
    K_z_half = sf_floatalloc(nzpad);
    for (i=0;i<nzpad;i++){
        K_z[i] = 1.0f;
        K_z_half[i] = 1.0f;
    }
    CPML_profile(d_x, K_x, alpha_prime_x, a_x, b_x, d_x_half, K_x_half, alpha_prime_x_half, a_x_half, b_x_half,
                 d_z, K_z, alpha_prime_z, a_z, b_z, d_z_half, K_z_half, alpha_prime_z_half, a_z_half, b_z_half,
                 nxpad, nzpad, nb, dx, dz, VPPML, dt, f0);
    
    
    /* =========================================Wavefield=========================================== */
    /*FD coefficients*/
    hc = sf_floatalloc(7);
    memset(hc,0,7*sizeof(float));
    gamma = coef(hc, maxerror, fdorder, MYID);
    /*check stability and dispersion*/
    vmax = velfwd[0][0];
    vmin = veladj[0][0];
    for(iz=0; iz<nz; iz++){
	for(ix=0; ix<nx; ix++){
	    vmax = SF_MAX(velfwd[iz][ix],vmax);
	    vmin = SF_MIN(veladj[iz][ix],vmin);
	}
    }
    if(vmin/(2*f0*hc[0]) < dx || vmin/(2*f0*hc[0]) < dz){
        fprintf(stderr,"dx,dz should be smaller than %f\n",vmin/(2*f0*hc[0]));
        MPI_Finalize(); 
        exit(1);
    }
    if(dx/(sqrt(2)*vmax*gamma) < dt){
        fprintf(stderr,"dt should be smaller than %f\n",dx/(sqrt(2)*vmax*gamma));
        MPI_Finalize(); 
        exit(1);
    }
    //src wfd
    sp =sf_floatalloc2(nxpad+2*fdo, nzpad+2*fdo); memset(sp [0],0,(nzpad+2*fdo)*(nxpad+2*fdo)*sizeof(float));
    svz=sf_floatalloc2(nxpad+2*fdo, nzpad+2*fdo); memset(svz[0],0,(nzpad+2*fdo)*(nxpad+2*fdo)*sizeof(float));
    svx=sf_floatalloc2(nxpad+2*fdo, nzpad+2*fdo); memset(svx[0],0,(nzpad+2*fdo)*(nxpad+2*fdo)*sizeof(float));
    if(!sf_getint("do_srcwfd", &do_srcwfd)) do_srcwfd=1;
    tmp_srcwfd_dir = sf_charalloc(512);
    if( NULL==(tmp_srcwfd_dir=sf_getstring("tmp_srcwfd_dir")) ){ 
        sf_warning("Should assign a tmp folder to store the source wavefield"); 
        MPI_Finalize(); exit(1);
    }
    if(do_srcwfd){
        if(!sf_getint("seismo_cal", &seismo_cal)) seismo_cal=0;
        if(seismo_cal){
            cal_dir_prefix = sf_charalloc(512);
            cal_dir_prefix = sf_getstring("cal_dir_prefix");
        }
    }
    //rec wfd
    gp =sf_floatalloc2(nxpad+2*fdo, nzpad+2*fdo); memset(gp [0],0,(nzpad+2*fdo)*(nxpad+2*fdo)*sizeof(float));
    gvz=sf_floatalloc2(nxpad+2*fdo, nzpad+2*fdo); memset(gvz[0],0,(nzpad+2*fdo)*(nxpad+2*fdo)*sizeof(float));
    gvx=sf_floatalloc2(nxpad+2*fdo, nzpad+2*fdo); memset(gvx[0],0,(nzpad+2*fdo)*(nxpad+2*fdo)*sizeof(float));
    //cpml
    psi1_dp_dx = sf_floatalloc2(nxpad, nzpad); memset(psi1_dp_dx[0], 0, nzpad*nxpad*sizeof(float));
    psi1_dp_dz = sf_floatalloc2(nxpad, nzpad); memset(psi1_dp_dz[0], 0, nzpad*nxpad*sizeof(float));
    psi1_dvx_dx = sf_floatalloc2(nxpad, nzpad); memset(psi1_dvx_dx[0], 0, nzpad*nxpad*sizeof(float));
    psi1_dvz_dz = sf_floatalloc2(nxpad, nzpad); memset(psi1_dvz_dz[0], 0, nzpad*nxpad*sizeof(float));
    psi2_dp_dx = sf_floatalloc2(nxpad, nzpad); memset(psi2_dp_dx[0], 0, nzpad*nxpad*sizeof(float));
    psi2_dp_dz = sf_floatalloc2(nxpad, nzpad); memset(psi2_dp_dz[0], 0, nzpad*nxpad*sizeof(float));
    psi2_dvx_dx = sf_floatalloc2(nxpad, nzpad); memset(psi2_dvx_dx[0], 0, nzpad*nxpad*sizeof(float));
    psi2_dvz_dz = sf_floatalloc2(nxpad, nzpad); memset(psi2_dvz_dz[0], 0, nzpad*nxpad*sizeof(float));
    //snapshot
    if(!sf_getint("snaptime", &snaptime)) snaptime=0; /* snap is # of time-step increment: snap=0 -> no snap; snap=1,2,3,... -> dt,2dt,3dt,... */
    if(snaptime){
        snapdir = sf_charalloc(512);
        snapdir = sf_getstring("snapdir"); /*folder of snapshot files*/
        if(!sf_getint("snaptype", &snaptype)) snaptype=1;
    }
    //data
    dat_dir_prefix = sf_charalloc(512);
    if( NULL==(dat_dir_prefix=sf_getstring("dat_dir_prefix")) ){ 
          printf("No external data directory & name-prefix inputs"); 
          MPI_Finalize(); exit(1);
    }
    
    /* =========================================Image=========================================== */
    den=sf_floatalloc2(nx,nz); memset(den[0], 0, nz*nx*sizeof(float));
    num=sf_floatalloc2(nx,nz); memset(num[0], 0, nz*nx*sizeof(float));
    img_1src=sf_floatalloc2(nx,nz); memset(img_1src[0], 0, nz*nx*sizeof(float));
    img_nsrc=sf_floatalloc2(nx,nz); memset(img_nsrc[0], 0, nz*nx*sizeof(float));
    img_out=sf_floatalloc2(nz,nx); memset(img_out[0], 0, nz*nx*sizeof(float));
    imgName = sf_charalloc(512);
    if( NULL==(imgName=sf_getstring("imgName")) ){ 
          fprintf(stderr,"No external imgName inputs");
          MPI_Finalize(); exit(1); 
    }
    if(!sf_getint("snapimg", &snapimg)) snapimg=0;//clduan, how many dt as an interval
    if(snapimg){
        //img_1src_it=sf_floatalloc2(nx,nz); memset(img_1src_it[0], 0, nz*nx*sizeof(float));
        snapimg_dir = sf_charalloc(512);
        snapimg_dir = sf_getstring("snapimg_dir");
    }
    if(!sf_getint("out_img_byshot", &out_img_byshot)) out_img_byshot=0;
    if(out_img_byshot){
        img_byshot_dir = sf_charalloc(512);
        img_byshot_dir = sf_getstring("img_byshot_dir");
    }
    
    /* ==========================Shot parallel========================== */
    if (MYID==0) sf_warning("Using %d processes\n",NP);
    if (MYID==0 && nsrc!=NP) { sf_warning("SRC_MPI: source num must be equal to Process num"); MPI_Finalize(); exit(1); }
    
    for (isrc=MYID; isrc<nsrc; isrc+=NP)
    {
	isx = srcpos[0][isrc]+fdo;
	isz = srcpos[1][isrc]+fdo;
	
	//receiver locations
	sprintf(receiver_file,"%srec_shot%d.txt",recv_dir,isrc+1);
	recpos = receiver(&ng, receiver_file, MYID, dx, dz, nx, nz, nb);
	
	//observed data
	dobs = sf_floatalloc2(nt,ng); memset(dobs[0], 0, ng*nt*sizeof(float));
	sprintf(datName,"%s_p.shot%d.rsf",dat_dir_prefix,isrc+1);
	if(MYID==0) sf_warning("Process 0 is reading: %s",datName);
        fdat = sf_input(datName);
        if(!sf_histint  (fdat, "n1", &c3)) {fprintf(stderr,"No n1= in dobs\n");MPI_Finalize();exit(1);}
        if(!sf_histint  (fdat, "n2", &c4)) {fprintf(stderr,"No n2= in dobs\n");MPI_Finalize();exit(1);}
        if(!sf_histfloat(fdat, "d1", &c5)) {fprintf(stderr,"No d1= in dobs\n");MPI_Finalize();exit(1);}
        if(c3!=nt || c4!=ng || c5!=dt){printf("dobs nt,dt,ng wrong");MPI_Finalize();exit(1);}
        sf_floatread(dobs[0], nt*ng, fdat);
	
	if(do_srcwfd){
	    if(seismo_cal){
	        dcal = sf_floatalloc2(nt,ng); memset(dcal[0], 0, ng*nt*sizeof(float));
	    }
	    for(it=0; it<nt; it++)
	    {
	        if(MYID==0 && (it+1)%pstp==0)  sf_warning("Source wavefield (forward) time-step %d of %d\n",it+1,nt);
	        sp[isz][isx] += wlt[it];
	        step_forward(nxpad, nzpad, dx, dz, dt, fdo, hc, inv_rho_half_x, inv_rho_half_z, rhopad, velfwd_pad,
                             K_x, a_x, b_x, K_z, a_z, b_z, K_x_half, a_x_half, b_x_half, K_z_half, a_z_half, b_z_half, 
                             psi1_dp_dx, psi1_dp_dz, psi1_dvx_dx, psi1_dvz_dz, sp, svx, svz);
	        //save to disk the source wavefield
	        bndr_rw(false, sp, tmp_srcwfd_dir, it, isrc, fdo, nb, nx, nz);
	    
	        //snapshot
	        if(snaptime && (it+1)%snaptime==0) outsnap(snaptype, sp, svx, svz, isrc, it, nxpad, nzpad, nb, fdo, snapdir, ".sw");
	    
	        //seismogram (used to compare with the observed data)
	        if(seismo_cal){
                    for (irec=0; irec<ng; irec++)
                        dcal[irec][it] = sp[recpos[1][irec]+fdo][recpos[0][irec]+fdo];
                }

	    }
	    if(seismo_cal){
	        sprintf(calName,"%s_p.bin.shot%d",cal_dir_prefix,isrc+1);
	        fdcal=fopen(calName,"w");
	        for(i=0;i<ng;i++)
	            for(j=0;j<nt;j++){
		        fwrite(&dcal[i][j],sizeof(float),1,fdcal); }
	        fclose(fdcal);
	    }
	    if(MYID==0)  sf_warning("FORWARD DONE!");
	}

//===================================================================	
	for(it=nt-1; it>-1; it--)
	{	
	    
            if(MYID==0 && (it+1)%pstp==0)  sf_warning("Receiver wavefield (adjoint) time-step %d of %d\n",it+1,nt);
	    
	    /*receiver wavefield*/
	    for(irec=0; irec<ng; irec++){
	        igx = recpos[0][irec]+fdo;
	        igz = recpos[1][irec]+fdo;
	        gp[igz][igx] += dobs[irec][it];
	    }
	    step_forward(nxpad, nzpad, dx, dz, dt, fdo, hc, inv_rho_half_x, inv_rho_half_z, rhopad, veladj_pad,
                         K_x, a_x, b_x, K_z, a_z, b_z, K_x_half, a_x_half, b_x_half, K_z_half, a_z_half, b_z_half, 
                         psi2_dp_dx, psi2_dp_dz, psi2_dvx_dx, psi2_dvz_dz, gp, gvx, gvz);

	    /*snapshot*/
	    if(snaptime && (it+1)%snaptime==0) outsnap(snaptype, gp, gvx, gvz, isrc, it, nxpad, nzpad, nb, fdo, snapdir, ".gw");
	    
	    /*reload the source wavefield*/
	    bndr_rw(true, sp, tmp_srcwfd_dir, it, isrc, fdo, nb, nx, nz);
	    
	    /*imaging condition*/
	    cross_correlation(num, den, sp, gp, nb, fdo, nx, nz);
	    
	    /*Allreduce snapshot of the image at it-step (clduan)*/
	    if(snapimg && (it+1)%snapimg==0)  outsnapimg(num, den, it, MYID, nx, nz, snapimg_dir);
            
	} //end time loop
	if(MYID==0)  sf_warning("ADJOINT DONE!");

	/*final-image = cross_correlation / illumination*/
	for(iz=0; iz<nz; iz++){
	    for(ix=0; ix<nx; ix++){
	        img_1src[iz][ix] = num[iz][ix]/(den[iz][ix]+SF_EPS);
	    }
	}
	//output shot-by-shot image
	if(out_img_byshot){
            sprintf(tmp_img1,"%simg_shot%d.bin",img_byshot_dir,isrc+1);
            fimg1=fopen(tmp_img1,"wb");
            for(i=0;i<nx;i++){
                for(j=0;j<nz;j++){
                    fwrite(&img_1src[j][i],sizeof(float),1,fimg1);
                }
            }
            fclose(fimg1);
        }
    
    } //end shot parallel
    
    /* (global) shot stacking */
    MPI_Allreduce(&img_1src[0][0], &img_nsrc[0][0], nx*nz, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    if(MYID==0) sf_warning("SHOT STACKING DONE!");
    
    /* output image */
    matrix_transp_2D(img_out, img_nsrc, nx, nz, 2);
    if (MYID==0){
	frtm = sf_output(imgName);
	sf_putint   (frtm, "n1", nz);
	sf_putfloat (frtm, "d1", dz);
	sf_putfloat (frtm, "o1", oz);
	sf_putstring(frtm, "label1", "Depth");
	sf_putstring(frtm, "unit1" , "m");
	sf_putint   (frtm, "n2", nx);
	sf_putfloat (frtm, "d2", dx);
	sf_putfloat (frtm, "o2", ox);
	sf_putstring(frtm, "label2", "Distance");
	sf_putstring(frtm, "unit2" , "m");
	sf_floatwrite(&img_out[0][0], nz*nx, frtm);
    }
    

    free(*srcpos);free(srcpos);
    free(wlt);
    if(src_func==3) free(wlt_name);
    if(outsrc==1 && MYID==0) free(outdir_src);
    free(*recpos);free(recpos);
    free(*vv1);free(vv1);
    free(*vv2);free(vv2);
    free(*rr);free(rr);
    free(*velfwd);free(velfwd);
    free(*veladj);free(veladj);
    free(*rho);free(rho);
    free(*velfwd_pad);free(velfwd_pad);
    free(*veladj_pad);free(veladj_pad);
    free(*rhopad);free(rhopad);
    free(*inv_rho_half_x);free(inv_rho_half_x);
    free(*inv_rho_half_z);free(inv_rho_half_z);
    free(d_x);
    free(K_x);
    free(alpha_prime_x);
    free(a_x);
    free(b_x);
    free(d_z);
    free(K_z);
    free(alpha_prime_z);
    free(a_z);
    free(b_z);
    free(d_x_half);
    free(K_x_half);
    free(alpha_prime_x_half);
    free(a_x_half);
    free(b_x_half);
    free(d_z_half);
    free(K_z_half);
    free(alpha_prime_z_half);
    free(a_z_half);
    free(b_z_half);
    free(hc);
    free(*sp);free(sp);
    free(*svx);free(svx);
    free(*svz);free(svz);
    free(*gp);free(gp);
    free(*gvx);free(gvx);
    free(*gvz);free(gvz);
    free(*psi1_dp_dx);free(psi1_dp_dx);
    free(*psi1_dp_dz);free(psi1_dp_dz);
    free(*psi1_dvx_dx);free(psi1_dvx_dx);
    free(*psi1_dvz_dz);free(psi1_dvz_dz);
    free(*psi2_dp_dx);free(psi2_dp_dx);
    free(*psi2_dp_dz);free(psi2_dp_dz);
    free(*psi2_dvx_dx);free(psi2_dvx_dx);
    free(*psi2_dvz_dz);free(psi2_dvz_dz);
    if(snaptime) free(snapdir);
    free(dat_dir_prefix);
    free(tmp_srcwfd_dir);
    if(do_srcwfd && seismo_cal){
        free(cal_dir_prefix);
        free(*dcal);free(dcal);
    }
    if(snapimg){
        //free(*img_1src_it);free(img_1src_it);
        free(snapimg_dir);
    }
    if(out_img_byshot){
        free(img_byshot_dir);
    }
    free(*den);free(den);
    free(*num);free(num);
    free(*img_1src);free(img_1src);
    free(*img_nsrc);free(img_nsrc);
    free(*img_out);free(img_out);
    free(imgName);
    free(*dobs);free(dobs);

    /* End MPI program */
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    exit(0);
    
}

