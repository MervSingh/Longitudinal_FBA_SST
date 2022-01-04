#!/bin/sh


SUBJLIST=`cat /home/mervyns/tg69_scratch/Mervyn/longit-fba-S3/subject.txt`

for SUBJ in $SUBJLIST

do

subid=`echo $SUBJ|awk '{print $1}' FS=","`
waveid=`echo $SUBJ|awk '{print $2}' FS=","`

sbatch --export sub=${subid},wave=${waveid} --job-name response_refs job_rf_total.sh

done
