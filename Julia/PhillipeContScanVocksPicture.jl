using DifferentialEquations
using DelimitedFiles
using Random
using Statistics
using Serialization

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

function CalcOrderMatrix(Ordermatrix, sol, Oszis)
    for k=1:Oszis,l=1:Oszis
        if abs((sol(10000)[k]-sol(3000)[k])/7000-(sol(10000)[l]-sol(3000)[l])/7000)<=0.001
            Ordermatrix[k,l] = 1.
        else
            Ordermatrix[k,l] = 0
        end
    end
end

function CalcFrequencies(freqs,sol,t)
    for i in eachindex(freqs)
        freqs[i]=(sol(t[2])[i]-sol(t[1])[i])/(t[2]-t[1])
    end
end

function ContScanthroughSigma(Initialstate, omegas;Sigmas=collect(0:0.01:10),α=0.0*pi,Betamult=-53,β=Betamult/100*pi)
    Order = zeros(length(omegas),length(Sigmas))
    Oszis = 50
    state = Initialstate
    t=(0.,11000)
    Frequencies = zeros(Oszis)
    params = CreateParams2(omegas,[collect(1:50),],α=α ,β=β ,σ=Sigmas[1])
    for i in 1:length(Sigmas)
        prob = ODEProblem(FullDerivatives2, state, t, params)
        sol = DifferentialEquations.solve(prob, progress = false, dtmax = 0.1, saveat=[10000,11000])
        CalcFrequencies(Frequencies, sol, (10000,11000))
        open(string("./",Betamult,"_",seed,"_",Sigmas[i],"_state.serialjl"),"w") do io
            serialize(io,sol(t[2]))
        end
        Order[1:end,i] .= Frequencies
        if i < length(Sigmas)
            params = CreateParams2(omegas,[collect(1:50),],α=α,β=β,σ=Sigmas[i+1])
            state = sol(t[2])
        end
    end
    return Order
end

seed = parse(Int64, ARGS[1])
Betamult = parse(Int64,ARGS[2])
state,omegas = SetupState(seed,50)
Order = ContScanthroughSigma(state,omegas,β = Betamult*pi/100)
open(string("./",Betamult,"_",seed,"_Frequencies.txt"),"w")  do io
    writedlm(io,Order)
end
