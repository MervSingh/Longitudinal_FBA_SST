#!/bin/sh


# set paths

data_dir='/scratch/tg69/dwi_preproc'
work_dir='/fs03/tg69/Mervyn/longit-fba'

# set variables
subid="sub-${sub}"
sessid="ses-${wave}"

cd $data_dir

#create rf folders
mkdir -p -v $work_dir/subjects_controls/${subid}_${sessid}
 

cp $data_dir/$subid/$sessid/response_csf_unbiased.txt $work_dir/subjects_controls/${subid}_${sessid}/response_csf_unbiased.txt
cp $data_dir/$subid/$sessid/response_gm_unbiased.txt $work_dir/subjects_controls/${subid}_${sessid}/response_gm_unbiased.txt
cp $data_dir/$subid/$sessid/response_wm_unbiased.txt $work_dir/subjects_controls/${subid}_${sessid}/response_wm_unbiased.txt
