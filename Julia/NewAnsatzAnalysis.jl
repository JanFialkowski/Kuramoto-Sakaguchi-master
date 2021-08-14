using DifferentialEquations
using Random
using Statistics
using DelimitedFiles
using Plots
using QuadGK
using Roots

function CreateContParams_(n...;α=0*pi,β=-0.53*pi,σ=3.5,ϵ=0.01,lowerlimit=-0.25,upperlimit=0.25)
    N=length(n)
    mean = zeros(N)
    var = zeros(N)
    lowerbound = lowerlimit
    for i in 1:N
        upperbound = lowerbound+n[i]*(upperlimit-lowerlimit)
        mean[i]=(lowerbound+upperbound)/2
        var[i]=1/12*(lowerbound-upperbound)^2
        lowerbound = upperbound
    end
    params = (N=N, Var = var, mean = mean, α=α, β=β, σ=σ, ϵ=ϵ,n=[n...])
end
function CreateContParams(n1;α=0*pi,β=-0.53*pi,σ=3.5,ϵ=0.01,lowerlimit=-0.25,upperlimit=0.25)
    N=2
    M1=(n1-1)*0.25
    M2 = n1*0.25
    V1 = n1^2/(3. *4*4)
    V2 = (1-n1)^2/(3. *4*4)
    params = (N=N, Var = [V1,V2], mean = [M1,M2], α=α, β=β, σ=σ, ϵ=ϵ,n=[n1,1-n1])
end
function CreateContParams(n...;α=0*pi,β=-0.53*pi,σ=3.5,ϵ=0.01,lowerlimit=-0.25,upperlimit=0.25)
    if sum(n)==1
        CreateContParams_(n...,α=α,β=β,σ=σ,ϵ=ϵ,lowerlimit=-0.25,upperlimit=0.25)
    elseif sum(n)<1
        CreateContParams_(n...,1-sum(n),α=α,β=β,σ=σ,ϵ=ϵ,lowerlimit=-0.25,upperlimit=0.25)
    else
        error("Sum of relative sizes greater than 1")
    end
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
function ContOrder(state,params::NamedTuple)
    R = zeros(params.N)
    Order = 0
    for i in 1:params.N
        R[i]=ContOrder(state,params,i)
    end
    for i in 1:params.N
        for j in i+1:params.N
            Order += 2*params.n[i]*R[i]*params.n[j]*R[j]*cos(state[params.N+i]-state[params.N+j])
        end
        Order+=params.n[i]*params.n[i]*R[i]*R[i]
    end
    Order
end
function ContOrder(state,orders,params::NamedTuple)
    R = orders
    Order = 0
    for i in 1:params.N
        for j in i+1:params.N
            Order += 2*params.n[i]*R[i]*params.n[j]*R[j]*cos(state[params.N+i]-state[params.N+j])
        end
        Order+=params.n[i]*params.n[i]*R[i]*R[i]
    end
    Order
end
function SteadyStateTheta(theta::Number, p, ci)
    R2 = ContOrder(2*theta,p.n[ci])
    1-p.σ/2/p.Var[ci]/theta*(cos(theta*p.n[ci]/2)-R2)*p.n[ci]*R2*sin(p.β)
end
function SteadyStateOrders(p)
    output1 = find_zero(x->SteadyStateTheta(x,p,1),1)
    output2 = find_zero(x->SteadyStateTheta(x,p,2),1)
    ContOrder(2*output1,p.n[1]),ContOrder(2*output2,p.n[2])
end
function OmegaPrime(Orders,p)
    phalf = p.mean[1]-p.mean[2]+p.σ*cos(p.β)/2*(p.n[1]*(1-Orders[1]^2)-p.n[2]*(1-Orders[2]^2))
    phalf /= -2
    q = -p.σ*p.ϵ*sin(p.β)/2
    if phalf*phalf-q<0
        return NaN,NaN
    else
        return -phalf + sqrt(phalf*phalf-q),-phalf-sqrt(phalf*phalf-q)
    end
end
ns = collect(0.5:0.05:0.95)
kaputts = zeros(length(ns))
sigmas = collect(0.5:0.01:4)
Orders = zeros((length(ns),length(sigmas)+1))
Orders[:,1].=ns
for i in 1:length(ns)
    for s in 1:length(sigmas)
        params = CreateContParams(ns[i],σ=sigmas[s])
        Order = SteadyStateOrders(params)
        if isnan(OmegaPrime(Order,params)[1])
            Orders[i,s+1]=NaN
        else
            Orders[i,s+1]=sqrt(ns[i]^2*Order[1]^2+(1-ns[i])^2*Orders[2]^2)
        end
    end
end
writedlm("./Data/NewestAnsatzTwoClusterOrders.txt",Orders)
plot(sigmas,transpose(Orders[:,2:end]),xlabel="Sigma",ylabel="Second Orderparameter")
