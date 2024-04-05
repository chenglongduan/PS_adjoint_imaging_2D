#!/bin/bash

nsrc=81
rec_dir='/scratch/scratch/cxd170430/AC_PSAdj_imaging_2D/ys_synth/receiver/'
in_bin_dir='/scratch/scratch/cxd170430/AC_PSAdj_imaging_2D/ys_synth/data/S_29km/stalta/'
out_rsf_dir='/scratch/scratch/cxd170430/AC_PSAdj_imaging_2D/ys_synth/data/S_29km/stalta/rsf/'


for (( i=1; i<=$nsrc; i++ ))
do
    nrec=$(< "${rec_dir}rec_shot${i}.txt" wc -l)
    
    echo n1=2001 n2=$nrec d1=0.004 d2=1 o1=0 o2=0 data_format="native_float" in="${in_bin_dir}obs.bin.shot${i}" label1="Time" label2="Distance" unit1="s" unit2="m" > ${out_rsf_dir}refl_p.shot${i}.rsf

done
