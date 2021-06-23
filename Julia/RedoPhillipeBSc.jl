using DifferentialEquations
using DelimitedFiles
using Random
using Statistics

function Derivativesphi2(dstate, system,params,t)
    dstate[1:50] .= 0
    sina = sin.(system[1:params.N] .+params.α)
    cosa = cos.(system[1:params.N] .+params.α)
    sin0 = sin.(system[1:params.N])
    cos0 = cos.(system[1:params.N])
    for i = 1:params.N
        for j = 1:params.N
            dstate[i] -= system[params.N*i+j] * (sina[i]*cos0[j]-cosa[i]*sin0[j])
        end
        dstate[i] *= params.σ / params.N
        dstate[i] += params.omegas[i]
    end
end

function Derivativeskappa2(dstate, system,params,t)
    dstate[51:end] .= 0
    sinb = sin.(system[1:params.N] .+params.β)
    cosb = cos.(system[1:params.N] .+params.β)
    sin0 = sin.(system[1:params.N])
    cos0 = cos.(system[1:params.N])
    for i in 1:params.N
        for j in 1:params.N
            dstate[params.N*i+j] = -params.ϵ*(system[params.N*i+j]+(sinb[i]*cos0[j]-cosb[i]*sin0[j]))
        end
    end
end

function FullDerivatives2(dstate,system,params,t)
    Derivativesphi2(dstate,system,params,t)
    Derivativeskappa2(dstate,system,params,t)
end

function CreateParams2(omegas, indices; α=0.3*pi, β=0.23*pi, σ=3.5, ϵ=0.01)
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
    params = (omegas = omegas, indices = indices, Σ = S, mean = m, α=α, β=β, σ=σ, N=length(omegas), ϵ=ϵ)
end

function SetupState(seed,Oszis)
    Random.seed!(seed)
    omegas = sort!(0.5*rand(Oszis)-0.25 .*ones(Oszis))
    state= 2*rand(Oszis+Oszis*Oszis)
    state[1:Oszis].*=pi
    state[Oszis+1:end] .-= 1.
    return state, omegas
end

function CalcOrderMatrix(Ordermatrix, sol, Oszis;t1=3000,t2=10000)
    for k=1:Oszis,l=1:Oszis
        if abs((sol(t2)[k]-sol(t1)[k])/(t2-t1)-(sol(t2)[l]-sol(t1)[l])/(t2-t1))<=0.001
            Ordermatrix[k,l] = 1.
        else
            Ordermatrix[k,l] = 0
        end
    end
end

function ContScanthroughSigma(Initialstate, omegas;Sigmas=collect(0:0.01:10),α=0.0*pi,Betamult=-53,β=Betamult/100*pi,filename="Ordermatrix")
    Order = zeros(length(Sigmas))
    Oszis = 50
    state = Initialstate
    t=(0.,10000)
    Ordermatrix = zeros((Oszis,Oszis))
    params = CreateParams2(omegas,[collect(1:50),],α=α ,β=β ,σ=Sigmas[1])
    for i in 1:length(Sigmas)
        prob = ODEProblem(FullDerivatives2, state, t, params)
        sol = solve(prob, progress = false, dtmax = 0.1, saveat=[3000,10000])
        CalcOrderMatrix(Ordermatrix, sol, Oszis)
        open(string("./",Betamult,"_",seed,"_",Sigmas[i],"_",filename,".txt"),"w") do io
            writedlm(io,Ordermatrix)
        end
        Order[i]=mean(Ordermatrix)
        if i < length(Sigmas)
            params = CreateParams2(omegas,[collect(1:50),],α=α,β=β,σ=Sigmas[i+1])
        end
        state = sol(10000)
    end
    return Order, state
end

#=
#doshitbackwards()
#bla = TestODE()
Sigmas = collect(0:1.0:4.)
Order = zeros(length(Sigmas))
Oszis = 50
Random.seed!(3)
omegas = sort!(0.5*rand(Oszis)-0.25 .*ones(Oszis))
state= rand(Oszis+Oszis*Oszis)
state[1:50].*=2*pi
state[51:end] .-= 0.
t=(0.,10000)
Ordermatrix = zeros((Oszis,Oszis))
params = CreateParams2(omegas,[collect(1:50),],α=0.3*pi,β=-0.53*pi,σ=3.5)
prob = ODEProblem(FullDerivatives2, state, t, params)
system = Kuramotosystem(zeros(Oszis),omegas = omegas, α=0.3*pi,β=-0.53*pi, σ=3.5)
system.state[1:end]=prob.u0[1:end]
#fullsol = CustomRK4(system, (0,10000), dt=0.001, save_indices = collect(1:Oszis+Oszis^2), save_times=[3000,10000])
#fullsol=readdlm("./Data/Phillipe/HighAccCustomSolve.txt")
=#
seed = parse(Int64, ARGS[1])
Betamult = parse(Int64,ARGS[2])
state,omegas = SetupState(seed,50)
Random.seed!(0)
omegas = sort!(0.5*rand(50)-0.25 .*ones(50))
Order,Finalstate = ContScanthroughSigma(state,omegas,β = Betamult*pi/100,Sigmas=collect(0:0.01:10),filename = "Upsweep")
open(string("./",Betamult,"_",seed,"_UpSweepOrders.txt"),"w")  do io
    writedlm(io,Order)
end
state = Finalstate
Order,Finalstate = ContScanthroughSigma(state,omegas,β = Betamult*pi/100,Sigmas=collect(10:-0.01:0),filename = "Downsweep")
open(string("./",Betamult,"_",seed,"_DownSweepOrders.txt"),"w")  do io
    writedlm(io,Order)
end
