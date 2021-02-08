using DifferentialEquations
using Plots
using DelimitedFiles
using Random
include("KuramotoCustomIntegrator.jl")

function Derivativesphi(dstate, system,params,t)
    dstate[1:params.N] .= 0
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

function Derivativeskappa(dstate, system,params,t)
    dstate[params.N+1:end] .= 0
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

function FullDerivatives(dstate,system,params,t)
    Derivativesphi(dstate,system,params,t)
    Derivativeskappa(dstate,system,params,t)
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
    params = (omegas = omegas, indices = indices, Σ = S, mean = m, α=α, β=β, σ=σ, N=Oszis, ϵ=ϵ)
end

Oszis = 50
Random.seed!(314159265)
omegas = sort!(0.5*rand(Oszis)-0.25 .*ones(Oszis))
params = CreateParams(omegas,[collect(1:50),],α=0.3*pi,β=-0.53*pi,σ=2.5)
state = [2*pi*rand(Oszis),rand(Oszis*Oszis)]
t=(0.,10000.)
prob = ODEProblem(Derivativesphi, state, t, params)
sol = solve(prob, progress = true)
