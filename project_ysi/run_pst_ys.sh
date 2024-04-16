#! /bin/bash


#SBATCH -J test
#SBATCH -N 3
#SBATCH -n 81
#SBATCH -p Optane
#SBATCH -t 06:00:00

##module load madagascar
##module load mvapich2
##module load fftw/3.3.8

export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MV2_ENABLE_AFFINITY=0


par_rtm=(
	fdorder=8     # only 8
	maxerror=0    # 0=Taylor; 1=0.1%Holberg; 2=0.02%Holberg
	nb=30	
	f0=5.0
	VPPML=4500.0 
	printstep=200    # print every ? steps	
	velocity_fwd="./model/ys_homo/vp.rsf"  #un-padded model, rsf file
	velocity_adj="./model/ys_homo/vs.rsf"  #un-padded model, rsf file
	density="./model/ys_homo/rho.rsf"    #un-padded model, rsf file
	dt=0.004
	tt=8.004
	source_file="./source/source_81_ys.txt" #real physical location [nsrc;sx,sz]
	amp=4.8e3
	src_func=2 #1=Ricker 2=Gaussian 3=customize
	src_wlt="./wavelet/wavelet.rsf"  #rsf file
	time_delay=0.0
	outsrc=1
	outdir_src="./wavelet/3/"    #bin file
	recv_dir="./receiver/ys/"  # [dir/]rec_shot?.txt
	dat_dir_prefix="./data/ys_sta_lta/T/ps/rsf/refl" #rsf file([dir/prefix]_p.shot%d.rsf)
	do_srcwfd=1
	tmp_srcwfd_dir="/petastore/lumleylab/cxd170430/tmp_srcwfd/"
	seismo_cal=0
	cal_dir_prefix="./outdata/1/10hz"  #bin file
	snaptime=0 #how many dt as an interval
	snaptype=1  #1=pressure 2=vx 3=vz
	snapdir="./snapshot/1/"    #directory to store snapshots, end with slash (/), bin file
	imgName="./result/vt_ps_dt0.rsf"  #rsf file
	snapimg=0  #how many dt as an interval
	snapimg_dir="./snapimg/1/"
	out_img_byshot=0
	img_byshot_dir="./result/img_shot/"
	datapath="/petastore/ganymede/home/cxd170430/mako_codes/ACOUSTIC2D_RTM_NEW/project_twolayer/tmp_outbin/"  #absolute directory to store output binary file
)


prun -v ../src/x_acrtm_ps ${par_rtm[@]}

