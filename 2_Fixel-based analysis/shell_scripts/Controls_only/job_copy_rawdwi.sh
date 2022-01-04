#!/bin/bash
#this script creates symbolic links for each subject's wmfod and mask images

module load mrtrix3tissue/5.2.8

data_dir='/projects/tg69/bids'
work_dir='/home/mervyns/tg69_scratch/Mervyn/longit-fba'

cd $work_dir

#set variables
subid="sub-${sub}"
sessid="ses-${wave}"

cp -v $data_dir/$subid/$sessid/dwi/${subid}_${sessid}_acq-mB2800_dwi*.{nii.gz,bvec,bval} $work_dir/subjects_controls/${subid}_${sessid}/${subid}_${sessid}_acq-mB2800_dwi*.{nii.gz,bvec,bval}



