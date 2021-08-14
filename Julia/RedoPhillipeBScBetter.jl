using DifferentialEquations
using DelimitedFiles
using Random
using Statistics
using QuadGK

function Derivativesphi2(dstate, system,params,t)
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

function Derivativeskappa2(dstate, system,params,t)
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
function KuramotoOrder(state)
    Order = 0. + 0.0im
    for phi in state
        Order += exp(1im*phi)
    end
    abs(Order/length(state))
end
function ContScanthroughSigma(Initialstate, omegas; Oszis = 100, Sigmas=collect(0:0.01:10),
    α=0.0*pi,Betamult=-53,β=Betamult/100*pi,filename="Ordermatrix",t=(0.,1000.),n1=50)

    Ordermatrix = zeros((Oszis,Oszis))
    resetscan = true
    Order = zeros(length(Sigmas))
    state = Initialstate
    t=(0.,10000)
    params = CreateParams2(omegas,[collect(1:Oszis),],α=α ,β=β ,σ=Sigmas[1])
    for i in 1:length(Sigmas)
        prob = ODEProblem(FullDerivatives2, state, (0.,10000.), params)
        sol = solve(prob, progress = false, dtmax = 0.1, saveat=[3000,10000])
        CalcOrderMatrix(Ordermatrix, sol, Oszis)
        Order[i]=mean(Ordermatrix)
        try
            if findfirst(isequal(0),Ordermatrix)[1]==n1+1 && findlast(isequal(0),Ordermatrix)[1]==n1
                state=sol(10000)
            end
        catch
            state=sol(10000)
        end
        if i < length(Sigmas)
            params = CreateParams2(omegas,[collect(1:Oszis),],α=α,β=β,σ=Sigmas[i+1])
        end
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
function Sweepdata(initstate,omegas;Betamult=-53,sigmas=collect(0:0.01:10))
    Order,Finalstate = ContScanthroughSigma(initstate,omegas,β = Betamult*pi/100,Sigmas=sigmas,filename = "Upsweep",Oszis=length(omegas))
    open(string("./",Betamult,"_",seed,"_UpSweepOrders.txt"),"w")  do io
        writedlm(io,Order)
    end
    Order,Finalstate = ContScanthroughSigma(Finalstate,omegas,β = Betamult*pi/100,Sigmas=sigmas,filename = "Downsweep",Oszis=length(omegas))
    open(string("./",Betamult,"_",seed,"_DownSweepOrders.txt"),"w")  do io
        writedlm(io,Order)
    end
end
function BifDiagramScharFullSystem(n1,omegas;sigmas=collect(0:0.01:10))
    sigmasreverse = sigmas[end:-1:1]
    for i in 1:length(n1)
        n = n1[i]
        initstate = zeros(length(omegas)+length(omegas)^2)
        initstate[n+1:end].=0.3*pi
        for i in 1:length(omegas)
            for j in 1:length(omegas)
                if (i<=n&&j<=n)||(i>n&&j>n)
                    initstate[length(omegas)*i+j]=1
                end
            end
        end
        Orders,Finalstate=ContScanthroughSigma(initstate,omegas,n1=n,Sigmas=sigmas,Oszis=length(omegas))
        writedlm(string("./FullSystemBifDiagramUpSweep_",n,".txt"),Orders)
        Orders,Finalstate=ContScanthroughSigma(Finalstate,omegas,n1=n,Sigmas=sigmasreverse,Oszis=length(omegas))
        writedlm(string("./FullSystemBifDiagramDownSweep_",n,".txt"),Orders)
    end
end
function TurnDataintoArray(Files,ns)
    Order=readdlm(files[1])
    Orders = zeros(length(ns),length(Order)+length(ns[1]))
    for i in 1:length(ns)
        Order = readdlm(files[i])
        Orders[i,1:length(ns[i])].=ns[i]...
        Orders[i,length(ns[i])+1:end]=Order
    end
    Orders
end
seed = parse(Int64, ARGS[1])
Betamult = parse(Int64,ARGS[2])
state,omegas = SetupState(seed,50)
Random.seed!(0)
omegas = LinRange(-0.25,0.25,50)
BifDiagramScharFullSystem([seed],omegas)
