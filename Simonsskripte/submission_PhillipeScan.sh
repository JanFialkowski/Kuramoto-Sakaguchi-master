#!/bin/bash
#loading parameters
. parameters_submission_PhillipeScan.sh
mkdir -p /$outdir/
#write logfile
echo "working directory on each computer --> $work" >> /$outdir/log_PhillipeScan.txt
echo "data will be saved to  --> /$outdir/" >> /$outdir/log_PhillipeScan.txt
#Submit jobs
qsub -mem 2 -speed 3 -m e -r SIGMA,0,$SIGMAstep PhillipeScan.sh
