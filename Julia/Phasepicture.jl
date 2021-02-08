using Plots
using DelimitedFiles
using Random
include("KuramotoCustomIntegrator.jl")

Oszis = 400
cutoff = 0.8
s = 3.5
t = (0.,3000.)
Random.seed!(314159265)
omegas = sort!(2*rand(Oszis)-1. .*ones(Oszis))
Phis = zeros(Oszis)
for i in round(Int, cutoff*length(Phis), RoundDown):length(Phis)
    Phis[i]+=0.3*pi
end
System = Kuramotosystem(Phis,omegas=omegas,α=0.15*pi, β=-0.53*pi, σ = s)
for i in round(Int, cutoff*length(Phis), RoundDown):length(Phis)
    for j in 1:round(Int, cutoff*length(Phis),RoundDown)-1
        System.state[Oszis+Oszis*(i-1)+j]=0
        System.state[Oszis+Oszis*(j-1)+i]=0
    end
end
System.state[Oszis+1:end].+= 0.01*rand.(Oszis*Oszis)
full_sol = CustomRK4(System, t, save_times=collect(2900:3000))
open(string("./Data/3_5Phases.txt"),"w")  do io
    writedlm(io,full_sol)
end
