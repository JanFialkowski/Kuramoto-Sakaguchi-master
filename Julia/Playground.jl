using LightGraphs
using DifferentialEquations
using Plots
using DelimitedFiles
using Serialization
using Statistics

n = collect(1:50)
sigma = collect(0.:0.01:10.)
avgorders = zeros(length(sigma))
for i in 1:length(sigma)
    Cliques = readdlm(string("./Data/ContScan15/",sigma[i],"_15cont_Cliques.csv"))
    avgorders[i]=size(Cliques)[1]
end
plot(sigma,avgorders,xlabel="Sigma",ylabel="Number of clusters")
savefig("EquiW_ContScanClusters")
