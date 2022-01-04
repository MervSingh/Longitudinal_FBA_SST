#!/bin/bash

#this script preprocesses diffusion data and prepares it for probabilistic ssmt-csd

#load packages and set paths
data_dir='/fs03/tg69/dwi_preproc'
work_dir='/fs03/tg69/Mervyn/longit-fba-S3'

module load mrtrix3tissue/5.2.8

#set variables
subid="sub-${sub}"
sessid="ses-${wave}"

#begin job
echo "Working on ${subid} ${sessid}"

cd $data_dir

ss3t_csd_beta1 $data_dir/$subid/$sessid/dwi_denoised_unringed_preproc_unbiased_upsampled.mif $work_dir/total_group_average_response_wm.txt $work_dir/subjects_total/${subid}_${sessid}/wmfod.mif $work_dir/total_group_average_response_gm.txt $work_dir/subjects_total/${subid}_${sessid}/gm.mif $work_dir/total_group_average_response_csf.txt $work_dir/subjects_total/${subid}_${sessid}/csf.mif -mask $data_dir/$subid/$sessid/mask_unbiased_combined_upsampled.mif

cd $work_dir/subjects_total

mtnormalise ${subid}_${sessid}/wmfod.mif ${subid}_${sessid}/wmfod_norm.mif ${subid}_${sessid}/gm.mif ${subid}_${sessid}/gm_norm.mif ${subid}_${sessid}/csf.mif ${subid}_${sessid}/csf_norm.mif -mask $data_dir/$subid/$sessid/mask_unbiased_combined_upsampled.mif

fod2dec ${subid}_${sessid}/wmfod_norm.mif ${subid}_${sessid}/decfod.mif -mask $data_dir/$subid/$sessid/mask_unbiased_combined_upsampled.mif
 
