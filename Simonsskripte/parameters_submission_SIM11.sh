#parameters SIM11_MSF_iH0, use Float for start/end and Int for steps
H0start=0.7
H0end=0.8
H0steps=2
DH0start=70.0
DH0end=80.0
DH0steps=2
#sigma, overall coupling strength
sigma=0.35
#directory where condor is working somewhere else
work=CONDOR_SIM11_MSF_iH0/sigma${sigma}
#directory where everything is stored in the end
outdir=home/simon30002/SIM11_MSF_iH0_SAVE/data_for_BER20b_sigma${sigma}
#filename julia
file=MSF_FHN_SIM11.jl
#directory of file
progdir=home/simon3000/SIM11_MSF_iH0
