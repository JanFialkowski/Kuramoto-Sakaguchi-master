#!/bin/bash
#load parameters
. /home/itsunderepsilon/Codeprojects/Master/parameters_submission_PhillipeScan.sh
#Calculate Sigma, because -r option on qsub is weird for floats
digits=3
sigma=$(python -c "print(round($SIGMAstart+($SIGMAend-$SIGMAstart)/($SIGMAstep)*$QVARSIGMA,$digits))")
#echo hostnames
echo "Sigma = $sigma at $(hostname)" >> /$outdir/hosts.txt
#clean working directories
rm -r /users/itsunderepsilon/$work/Sigma_$sigma/
#create new working directory (one for each iteration)
mkdir -p /users/itsunderepsilon/$work/Sigma_$sigma/
#go to directory where julia file is located
cd /$progdir/
#copy julia file
cp $file /users/itsunderepsilon/$work/Sigma_$sigma/
#go to working direeeeectory 
cd /users/itsunderepsilon/$work/Sigma_$sigma/
#execute julia file
julia $file $sigma
#copy results to my computer
cp /users/itsunderepsilon/$work/Sigma_$sigma/* /$outdir/
#log successfull copying
echo "data sigma = $sigma written to /$outdir/" >> /$outdir/log_PhillipeScan.txt
#clean working directory
rm -r /users/itsunderepsilon/$work/Sigma_$sigma/
