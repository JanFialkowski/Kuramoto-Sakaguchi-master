#!/bin/bash
#load parameters
. /home/itsunderepsilon/Codeprojects/Master/parameters_submission_PhillipeScan.sh
#echo hostnames
echo "at $(hostname)" >> /$outdir/hosts.txt
#clean working directories
rm -r /users/itsunderepsilon/$work/Sigma/
#create new working directory (one for each iteration)
mkdir -p /users/itsunderepsilon/$work/Sigma/
#go to directory where julia file is located
cd /$progdir/
#copy julia file
cp $file /users/itsunderepsilon/$work/Sigma/
#go to working direeeeectory 
cd /users/itsunderepsilon/$work/Sigma/
#execute julia file
julia $file
#copy results to my computer
cp /users/itsunderepsilon/$work/Sigma/* /$outdir/
#log successfull copying
echo "data written to /$outdir/" >> /$outdir/log_PhillipeScan.txt
#clean working directory
rm -r /users/itsunderepsilon/$work/Sigma/
