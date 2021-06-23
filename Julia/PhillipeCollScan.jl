using DifferentialEquations
using Random
using Statistics
using DelimitedFiles
using QuadGK

function CreateContParams(n1;α=0*pi,β=-0.53*pi,σ=3.5,ϵ=0.01)
    M1=(n1-1)*0.25
    M2 = n1*0.25
    V1 = n1^2/(3. *4*4)
    V2 = (1-n1)^2/(3. *4*4)
    params = (N=2, Var = [V1,V2], mean = [M1,M2], α=α, β=β, σ=σ, ϵ=ϵ,n=[n1,1-n1])
end
function CreateContParams(n1,n2;α=0*pi,β=-0.53*pi,σ=3.5,ϵ=0.01)
    if n1+n2!=1
        error("n1+n2!=1")
    end
    CreateContparams(n1,α=α,β=β,σ=σ,ϵ=ϵ)
end
function ContOrder(state,p,ci)
    ContOrder(state[ci],p.n[ci])
end
function ContOrder(theta,n)
    if theta == 0
        1
    else
        4/(theta*n)*sin(theta*n/4)
    end
end
function ContDerivs(du,u,p::NamedTuple,t)
    R1 = ContOrder(u,p,1)
    R2 = ContOrder(u,p,2)
    du[1]=1+p.σ/p.Var[1]/u[1]*(cos(u[1]*p.n[1]/4)-R1)*(u[2*p.N+1]*p.n[1]*R1*cos(p.α)+u[2*p.N+2]*p.n[2]*R2*cos(u[p.N+1]-u[p.N+2]+p.α))
    du[2]=1+p.σ/p.Var[2]/u[2]*(cos(u[2]*p.n[2]/4)-R2)*(u[2*p.N+4]*p.n[2]*R2*cos(p.α)+u[2*p.N+3]*p.n[1]*R1*cos(u[p.N+2]-u[p.N+1]+p.α))
    du[3]=p.mean[1]-p.σ*p.n[1]*u[2*p.N+1]*sin(p.α)*R1^2-p.σ*p.n[2]*u[p.N*2+2]*R1*R2*sin(u[p.N+1]-u[p.N+2]+p.α)
    du[4]=p.mean[2]-p.σ*p.n[2]*u[2*p.N+4]*sin(p.α)*R2^2-p.σ*p.n[1]*u[p.N*2+3]*R1*R2*sin(u[p.N+2]-u[p.N+1]+p.α)
    du[5]=-p.ϵ*(u[p.N*2+1]+R1*R1*sin(p.β))
    du[8]=-p.ϵ*(u[p.N*2+4]+R2*R2*sin(p.β))
    du[6]=-p.ϵ*(u[p.N*2+2]+R1*R2*sin(u[p.N+1]-u[p.N+2]+p.β))
    du[7]=-p.ϵ*(u[p.N*2+3]+R1*R2*sin(u[p.N+2]-u[p.N+1]+p.β))
    if u[1]==0
        du[1]=1.
    end
    if u[2]==0
        du[2]=1.
    end
    du
end
function ContOrder(state,params::NamedTuple)
    R1 = ContOrder(state,params,1)
    R2 = ContOrder(state,params,2)
    n1 = params.n[1]
    n2 = params.n[2]
    n1^2*R1^2+n2^2*R2^2+n1*n2*R1*R2*cos(state[3]-state[4])
end
function ContScanthroughSigma(Initialstate, n1;Sigmas=collect(0:0.01:10),α=0.0*pi,Betamult=-53,β=Betamult/100*pi,filename="Ordermatrix")
    Order = zeros(length(Sigmas))
    state = Initialstate
    t=(0.,10000)
    params = CreateContParams(n1,α=α ,β=β ,σ=Sigmas[1])
    for i in 1:length(Sigmas)
        prob = ODEProblem(ContDerivs, state, (0.,5000), params)
        sol = solve(prob, progress = false, saveat=5000,dtmax=0.01)
        prob = ODEProblem(ContDerivs,sol(5000),t,params)
        sol = solve(prob)
        Order[i]=quadgk(x->ContOrder(sol(x),params),0,10000)[1]/10000.
        if i < length(Sigmas)
            params = CreateContParams(n1,α=α,β=β,σ=Sigmas[i+1])
        end
        if rem(Sigmas[i],1)==0
            writedlm(string("./",n1*50,"_",seed,"_",Sigmas[i],"_",filename,".txt"),sol)
        end
        state = sol(10000).+0.01.*rand(8)
    end
    return Order, state
end

seed = 138513#parse(Int64, ARGS[1])
n1mult = 27#parse(Int64,ARGS[2])
for n1mult in n1mult
    Random.seed!(seed)
    state = rand(8)
    state[3:4].*=2*pi
    n1 = n1mult/50.
    Order,Finalstate = ContScanthroughSigma(state, n1, β = -53*pi/100,Sigmas=collect(0.5:0.01:10),filename = "Upsweep")
    open(string("./",n1mult,"_",seed,"_UpSweepOrders.txt"),"w")  do io
        writedlm(io,Order)
    end
    state = Finalstate
    Order,Finalstate = ContScanthroughSigma(state,n1,β = -53*pi/100,Sigmas=collect(10:-0.01:0),filename = "Downsweep")
    open(string("./",n1mult,"_",seed,"_DownSweepOrders.txt"),"w")  do io
        writedlm(io,Order)
    end
end
