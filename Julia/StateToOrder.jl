using DifferentialEquations
using DelimitedFiles
using Random
using Statistics
using Serialization
using QuadGK

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

function order(Phis::Array{Float64,1})
    r = 0.0 + 0im
    for i = 1:length(Phis)
        r += exp(1im * Phis[i])
    end
    r /= length(Phis)
    r = abs(r)
end

function TimeAvgOrder(sol,t1,t2,n)
    Order,Error = quadgk(t->order(2 .*sol(t)[1:n]),t1,t2)
    Order/(t2-t1),Error
end

function RunState(state,params,t)
    prob = ODEProblem(FullDerivatives2,state,t,params)
    sol = solve(prob,dtmax=0.1,save_idxs=collect(1:length(params.omegas)))
end

function CalculateOrder(state,omegas,sigma;α=0,β=-0.53*pi)
    params = CreateParams2(omegas,[collect(1:length(omegas)),],α=α,β=β,σ=sigma)
    t = (0.,1000.)
    sol = RunState(state,params,t)
    Order = TimeAvgOrder(sol,t[1],t[2],length(omegas))[1]
end

function CreateOmegas(seed)
    state,omegas = SetupState(seed,50)
    omegas
end

function ReadState(path,sigma,seed)
    state = Array{Float64}(undef,1)
    open(string(path,"/-53_$(seed)_$(sigma)_state.serialjl")) do io
        state = deserialize(io)
    end
    state
end

function OrderFromPath(path,sigma,seed,t=(0.,1000.);oseed=seed)
    state = ReadState(path,sigma,seed)
    omegas = CreateOmegas(oseed)
    Order = CalculateOrder(state,omegas,sigma)
end

function ThatThing(path,seed,sigmas = collect(0:0.01:10);oseed=seed)
    Orders = zeros(length(sigmas))
    for i in 1:length(sigmas)
        Orders[i] = OrderFromPath(path,sigmas[i],seed,oseed=oseed)
    end
    Orders
end

function SetupState(seed,Oszis)
    Random.seed!(seed)
    omegas = sort!(0.5*rand(Oszis)-0.25 .*ones(Oszis))
    state= 2*rand(Oszis+Oszis*Oszis)
    state[1:Oszis].*=pi
    state[Oszis+1:end] .-= 1.
    return state, omegas
end

function Orderparameter(frequencies)
    Matrix = zeros(50,50)
    for i in 1:50
        for j in 1:50
            Matrix[i,j]=abs(frequencies[i]-frequencies[j])<0.001
        end
    end
    Matrix
end
function RepairFrequencies()
    Synchs = readdlm("./Data/Eckdaten/Synchronizationparameter_101runs_50Oszis.txt")
    writedlm("./Data/Eckdaten/Synchronizationparameter_101runs_50OszisOLD.txt",Synchs)
    sigmas = collect(0:0.01:10)
    bla = [0]
    for i in 0:100
        index = findlast(x->x<1,Synchs[i+1,1:end])
        if Synchs[index-1]>Synchs[index]
            bla = [bla index]
            state = ReadState("./Data/50OszisFrequencies_randomseed/Data",sigmas[index],i)
            omegas = CreateOmegas(i)
            params = CreateParams2(omegas,[collect(1:50),],α=0,β=-0.53*pi,σ=sigmas[index])
            sol = RunState(state,params,(0.,10000.))
            newfreqs = (sol(10000)[1:end] .- sol(9000)[1:end]) ./ 1000
            freqs = readdlm("./Data/50OszisFrequencies_randomseed/Data/-53_$(i)_Frequencies.txt")
            writedlm("./Data/50OszisFrequencies_randomseed/Data/-53_$(i)_FrequenciesOLD.txt",freqs)
            freqs[1:end,index] .= newfreqs
            writedlm("./Data/50OszisFrequencies_randomseed/Data/-53_$(i)_Frequencies.txt",freqs)
            Synchs[i+1,index]=mean(Orderparameter(newfreqs))
        end
    end
    writedlm("./Data/Eckdaten/Synchronizationparameter_101runs_50Oszis.txt",Synchs)
end

path = ARGS[1]
seed = parse(Int64,ARGS[2])
omegaseed=seed
if length(ARGS)==3
    omegaseed = parse(Int64,ARGS[3])
end
Orders = ThatThing(path,seed,oseed=omegaseed)
open(string("./-53_",seed,"_KuramotoOrders.txt"),"w")  do io
    writedlm(io,Orders)
end
