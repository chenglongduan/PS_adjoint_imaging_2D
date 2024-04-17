#!/bin/bash


# ------------concat all the SU data-----------------

#dir_su=/data/Chenglong/Yellowstone_vibroseis/su_data_sta_lta
#dir_comb=/data/Chenglong/Yellowstone_vibroseis/su_data_sta_lta
#dir_su=/data/Chenglong/Yellowstone_vibroseis/su_data_envelope
#dir_comb=/data/Chenglong/Yellowstone_vibroseis/su_data_envelope
dir_su=/data/Chenglong/Yellowstone_vibroseis/su_data_vel_traces
dir_comb=/data/Chenglong/Yellowstone_vibroseis/su_data_vel_traces



cat $dir_su/0911/*.R.su $dir_su/0912/*.R.su $dir_su/0913/*.R.su $dir_su/0914/*.R.su > $dir_comb/combine.R.su
cat $dir_su/0911/*.T.su $dir_su/0912/*.T.su $dir_su/0913/*.T.su $dir_su/0914/*.T.su > $dir_comb/combine.T.su
cat $dir_su/0911/*.Z.su $dir_su/0912/*.Z.su $dir_su/0913/*.Z.su $dir_su/0914/*.Z.su > $dir_comb/combine.Z.su
cat $dir_su/0911/*.RT.su $dir_su/0912/*.RT.su $dir_su/0913/*.RT.su $dir_su/0914/*.RT.su > $dir_comb/combine.RT.su
