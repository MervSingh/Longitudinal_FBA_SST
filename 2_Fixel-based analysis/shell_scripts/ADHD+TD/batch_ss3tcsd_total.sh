#!/bin/sh


SUBJLIST=`cat /home/mervyns/tg69_scratch/Mervyn/longit-fba-S3/subject.txt`

for SUBJ in $SUBJLIST

do

subid=`echo $SUBJ|awk '{print $1}' FS=","`
waveid=`echo $SUBJ|awk '{print $2}' FS=","`

sbatch --export sub=${subid},wave=${waveid} --job-name job_ss3tcsd_total --mem=4G --account=tg69 --time=10:00:00 job_script_ss3tcsd_total.sh

done
