#!/bin/sh


SUBJLIST=`cat /home/mervyns/tg69_scratch/Mervyn/longit-fba-S3/subject.txt`

for SUBJ in $SUBJLIST

do

subid=`echo $SUBJ|awk '{print $1}' FS=","`
waveid=`echo $SUBJ|awk '{print $2}' FS=","`

sbatch --export sub=${subid},wave=${waveid} --job-name copymasks --account=tg69 --time=2:00:00 job_copy_masks_total.sh

done
