#!/bin/bash
#load parameters
. /home/simon3000/SIM11_MSF_iH0/parameters_submission_SIM11.sh
#parameter transform
digits=3 #accuracy paramter calculation (~stepsize in parameters_submission_SIM11.sh)
H0=$(python -c "print(round($H0start+($H0end-$H0start)/($H0steps-1)*$QVARH0,$digits))")
DH0=$(python -c "print(round($DH0start+($DH0end-$DH0start)/($DH0steps-1)*$QVARDH0,$digits))")
#echo hostnames
echo "H0 = ${H0}, DH0 = ${DH0} at $(hostname)" >> /$outdir/hosts.txt
#clean working directories
rm -r /users/simon3000/$work/H0${H0}_DH0${DH0}/
#create new working directory (one for each iteration)
mkdir -p /users/simon3000/$work/H0${H0}_DH0${DH0}/
#go to directory where julia file is located
cd /$progdir/
#copy julia file
cp $file /users/simon3000/$work/H0${H0}_DH0${DH0}/
#go to working direeeeectory 
cd /users/simon3000/$work/H0${H0}_DH0${DH0}/
#execute julia file
julia $file $sigma $H0 $DH0
#copy results to my computer
cp /users/simon3000/$work/H0${H0}_DH0${DH0}/* /$outdir/
#log successfull copying
echo "data H0 = $H0 DH0 = $DH0 written to /$outdir/" >> /$outdir/log_submission_SIM11_MSF_iH0.txt
#clean working directory
rm -r /users/simon3000/$work/H0${H0}_DH0${DH0}/
