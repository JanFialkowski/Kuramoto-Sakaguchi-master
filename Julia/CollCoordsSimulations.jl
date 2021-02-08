include("CollCoordsFunctions.jl")
using DifferentialEquations
using Random
using DelimitedFiles

save_times=collect(2500:0.25:3000)
solution = zeros
Sigmas = collect(2.5:0.1:4.0)
Oszis=400
cutoff = 0.8
for s in Sigmas
    Random.seed!(314159265)
    omegas = sort!(2*rand(Oszis)-1. .*ones(Oszis))
    indices = [collect(1:round(Int, cutoff*length(omegas), RoundDown)-1),collect(round(Int, cutoff*length(omegas), RoundDown):Oszis)]
    params = CreateParams(omegas,indices, α=0.15*pi, β=-0.53*pi, σ = s)
    u = [0.0, 0.0, 0., 0.3*pi, 1., 0., 0., 1.]
    t = (0.,3000.)
    prob = ODEProblem(derivativesCollCoords, u, t, params)
    sol = solve(prob, progress = true)
    solution = zeros(length(save_times),8)
    for i in 1:length(save_times)
        solution[i,1:end] = sol(save_times[i])
    end
    open(string("./Julia/CollSims/",s,".txt"),"w")  do io
        writedlm(io,solution)
    end
end
