#!/bin/bash
#this script creates symbolic links for each subject's wmfod and mask images

module load mrtrix3tissue/5.2.8

data_dir='/fs03/tg69/dwi_preproc'
work_dir='/home/mervyns/tg69_scratch/Mervyn/longit-fba-S3'

cd $work_dir

#set variables
subid="sub-${sub}"
sessid="ses-${wave}"

cp $data_dir/$subid/$sessid/mask_unbiased_combined_upsampled.mif $work_dir/subjects_total/${subid}_${sessid}/mask_unbiased_combined_upsampled.mif

