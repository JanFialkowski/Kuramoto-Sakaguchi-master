#!/bin/bash
#loading parameters
. ./parameters_submission_SIM11.sh
mkdir -p /$outdir/
#write logfile
echo "working directory on each computer --> $work" >> /$outdir/log_submission_SIM11_MSF_iH0.txt
echo "data will be saved to  --> /$outdir/" >> /$outdir/log_submission_SIM11_MSF_iH0.txt
#Submit jobs
qsub -mem 2 -speed 3 -m n -r H0,0,$((H0steps-1)) -r DH0,0,$((DH0steps-1)) SIM11_FHN_iH0.sh
