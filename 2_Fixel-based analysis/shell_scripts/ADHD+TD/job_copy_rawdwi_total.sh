#!/bin/bash
#this script creates symbolic links for each subject's wmfod and mask images

module load mrtrix3tissue/5.2.8

copy_dir='/home/mervyns/tg69/bids'
work_dir='/home/mervyns/tg69_scratch/Mervyn/longit-fba-S3/subjects_total'

cd $copy_dir

#set variables
subid="sub-${sub}"
sessid="ses-${wave}"

cp -v $subid/$sessid/dwi/"*acq-mB2800*_dwi.nii.gz" $work_dir/${subid}_${sessid}/"*acq-mB2800*_dwi.nii.gz"

cp -v $subid/$sessid/dwi/"*acq-mB2800*_dwi.bvec" $work_dir/${subid}_${sessid}/"*acq-mB2800*_dwi.bvec"

cp -v $copy_dir/$subid/$sessid/dwi/"*acq-mB2800*_dwi.bval" $work_dir/${subid}_${sessid}/"*acq-mB2800*_dwi.bval"  

