using Random
using DifferentialEquations
using DelimitedFiles
using Statistics

function createsins(u, params; func = sin, extra = 0)
    sins = zeros(length(params.omegas))
    for i in 1:params.N
        for k in params.indices[i]
            sins[k] = func(u[i]*(params.omegas[k]-params.mean[i])+u[params.N+i]+extra)
        end
    end
    return sins
end

function derivativestheta(du, u,  params, t; sina=createsins(u,params, extra = params.α), cosa = createsins(u,params, func = cos, extra = params.α), sin0 = createsins(u,params), cos0 = createsins(u,params, func = cos))
    du[1:params.N].=0.0
    for i in 1:params.N#Index des clusters von theta
        for j in 1:params.N#Index der anderen Cluster
            for k in params.indices[i]#Oszillatorindex im Thetacluster
                for l in params.indices[j]#Oszillatorindex im anderen Cluster
                    du[i]-=u[2*params.N+params.N*(i-1)+j]*(params.omegas[k]-params.mean[i])*(sina[k]*cos0[l]-cosa[k]*sin0[l])
                end
            end
        end
        du[i]*=params.σ/length(params.omegas)/params.Σ[i]/size(params.indices[i])[1]
    end
    du[1:params.N].+=1.0
end

function derivativesf(du, u,  params, t; sina=createsins(u,params, extra = params.α), cosa = createsins(u,params, func = cos, extra = params.α), sin0 = createsins(u,params), cos0 = createsins(u,params, func = cos))
    du[params.N+1:2*params.N].=0.0
    for i in 1:params.N#Index des clusters von F
        for j in 1:params.N#Index der anderen Cluster
            for k in params.indices[i]#Oszillatorindex im F-Cluster
                for l in params.indices[j]#Oszillatorindex im anderen Cluster
                    du[params.N+i]-=u[2*params.N+params.N*(i-1)+j]*(sina[k]*cos0[l]-cosa[k]*sin0[l])
                end
            end
        end
        du[params.N+i]*=params.σ/length(params.omegas)/size(params.indices[i])[1]
    end
    du[params.N+1:2*params.N].+=params.mean
end

function derivativeskappa(du, u,  params, t; sinb=createsins(u,params, extra = params.β), cosb = createsins(u,params, func = cos, extra = params.β), sin0 = createsins(u,params), cos0 = createsins(u,params, func = cos))
    du[2*params.N+1:end].=0.0
    for i in 1:params.N#Index des ersten clusters
        for j in 1:params.N#Index des zweiten Cluster
            for k in params.indices[i]#Oszillatorindex im ersten Clsuter
                for l in params.indices[j]#Oszillatorindex im anderen Cluster
                    du[2*params.N+params.N*(i-1)+j]+=(sinb[k]*cos0[l]-cosb[k]*sin0[l])
                end
            end
            du[2*params.N+params.N*(i-1)+j]/=(size(params.indices[i])[1]*size(params.indices[j])[1])
        end
    end
    du[2*params.N+1:end].+=u[2*params.N+1:end]
    du[2*params.N+1:end].*=-params.ϵ
end

function derivativesCollCoords(du, u,  params, t)
    sina = createsins(u,params, extra = params.α)
    cosa = createsins(u,params, func = cos, extra = params.α)
    sin0 = createsins(u,params)
    cos0 = createsins(u,params, func = cos)
    derivativestheta(du, u, params, t, sina=sina, sin0=sin0, cosa=cosa, cos0=cos0)
    derivativesf(du, u, params, t, sina=sina, sin0=sin0, cosa=cosa, cos0=cos0)
    derivativeskappa(du, u, params, t, sin0=sin0, cos0=cos0)
end

function ThetastoPhis(u,params)
    Phis = zeros(size(params.omegas))
    for i in eachindex(params.indices)
        for j in params.indices[i]
            Phis[j]=u[i]*(params.omegas[j]-params.mean[i])+u[params.N+i]
        end
    end
    return Phis
end

function CreateParams(omegas, indices; α=0.3*pi, β=0.23*pi, σ=3.5, ϵ=0.01)
    n = size(indices)[1]
    m = zeros(n)
    S = zeros(n)
    for i in 1:n
        for j in indices[i]
            m[i]+=omegas[j]
        end
        m[i]/=size(indices[i])[1]
    end
    for i in 1:n
        for j in indices[i]
            S[i]+=(omegas[j]-m[i])^2
        end
        S[i]/=size(indices[i])[1]
    end
    params = (omegas = omegas, indices = indices, Σ = S, mean = m, α=α, β=β, σ=σ, N=n, ϵ=ϵ)
end

function CalcOrderMatrix(Ordermatrix, sol, params::NamedTuple)
    Phi10000=ThetastoPhis(sol(10000),params)
    Phi3000=ThetastoPhis(sol(3000),params)
    for k=1:length(params.omegas),l=1:length(params.omegas)
        if abs((Phi10000[k]-Phi3000[k])/7000-(Phi10000[l]-Phi3000[l])/7000)<=0.001
            Ordermatrix[k,l] = 1.
        else
            Ordermatrix[k,l] = 0
        end
    end
end
function ContScanthroughSigma(Initialstate, omegas;Sigmas=collect(0:0.01:10),α=0.0*pi,Betamult=-53,β=Betamult/100*pi,filename="Ordermatrix")
    Order = zeros(length(Sigmas))
    Oszis = length(omegas)
    state = Initialstate
    t=(0.,10000)
    Ordermatrix = zeros((Oszis,Oszis))
    params = CreateParams(omegas,[collect(1:27),collect(28:50)],α=α ,β=β ,σ=Sigmas[1])
    for i in 1:length(Sigmas)
        prob = ODEProblem(derivativesCollCoords, state, t, params)
        sol = solve(prob, progress = false, dtmax = 0.1, saveat=[3000,10000])
        CalcOrderMatrix(Ordermatrix, sol, params)
        open(string("./",Betamult,"_",seed,"_",Sigmas[i],"_",filename,".txt"),"w") do io
            writedlm(io,Ordermatrix)
        end
        Order[i]=mean(Ordermatrix)
        if i < length(Sigmas)
            params = CreateParams(omegas,[collect(1:27),collect(28:50)],α=α,β=β,σ=Sigmas[i+1])
        end
        state = sol(10000)
    end
    return Order, state
end
#=
seed = 0#parse(Int64, ARGS[1])
Betamult = -53#parse(Int64,ARGS[2])
Random.seed!(seed)
omegas=sort(0.5.*rand(50)-0.25.*ones(50))
indices=[collect(1:27),collect(28:50)]
state = rand(8)
state[3:4] .*= pi*2
Order,Finalstate = ContScanthroughSigma(state,omegas,β = Betamult*pi/100,Sigmas=collect(0:0.01:1),filename = "Upsweep")
open(string("./",Betamult,"_",seed,"_UpSweepOrders.txt"),"w")  do io
    writedlm(io,Order)
end
state = Finalstate
Order,Finalstate = ContScanthroughSigma(state,omegas,β = Betamult*pi/100,Sigmas=collect(1:-0.01:0),filename = "Downsweep")
open(string("./",Betamult,"_",seed,"_DownSweepOrders.txt"),"w")  do io
    writedlm(io,Order)
end
=#
t=(0.,10000.)
state = rand(8)
state[3:4] .*=pi*2
Random.seed!(0)
omegas = sort(0.5.*rand(50)-0.25.*ones(50))
indices=[collect(1:27),collect(28:50)]
params = CreateParams(omegas,indices, α=0., β=-0.53*pi,σ=10.)
Problem = ODEProblem(derivativesCollCoords, state,t,params)
sol = DifferentialEquations.solve(Problem,dtmax=0.1,saveat=[3000,10000])
open("./CollSim_0_28.txt","w") do io
    #writedlm(io,sol.u)
end
