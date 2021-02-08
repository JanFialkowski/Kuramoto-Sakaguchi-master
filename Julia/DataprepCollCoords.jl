using DelimitedFiles
using Statistics
using Random
include("CollCoordsFunctions.jl")
include("KuramotoCustomIntegrator.jl")

Oszis=400
Sigmas = collect(2.5:0.1:4.0)
indices = (collect(1:319),collect(320:400))
Random.seed!(314159265)
omegas = sort!(2*rand(Oszis)-1. .*ones(Oszis))
for s in Sigmas
    params = CreateParams(omegas,indices, α=0.15*pi, β=-0.53*pi, σ = s)
    Data = readdlm(string("./Julia/CollSims/",s,".txt"),Float64)
    Omegas = [0.,0.]
    Omegas[1] = (Data[end,3]-Data[end-401,3])/100
    Omegas[2] = (Data[end,4]-Data[end-401,4])/100
    Kappas = Data[end-400:end,5:8]
    Orders = [zeros(401),zeros(401)]
    for t in 1:401
        Phis = ThetastoPhis(Data[t,:],params)
        Orders[1][t] = order(Phis[1:319])
        Orders[2][t] = order(Phis[320:400])
    end
    open(string("./Data/CollOmegas_",s,".txt"),"w") do io
        writedlm(io,Omegas)
    end
    open(string("./Data/CollOrders_",s,".txt"),"w") do io
        writedlm(io,Orders)
    end
    open(string("./Data/CollKappas_",s,".txt"),"w") do io
        writedlm(io,Kappas)
    end
end
