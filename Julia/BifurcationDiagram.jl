using DifferentialEquations, Plots, ColorSchemes
using Random
using Statistics
using DelimitedFiles
using QuadGK

function SteadyClusterStateCondition(state,p,ci)
    SteadyClusterStateCondition(state[ci],p,ci)
end
function SteadyClusterStateCondition(theta::Number,p,ci)
    R = ContOrder(theta,p.n[ci])
    1-p.σ/theta/p.Var[ci]*(cos(theta*p.n[ci]/4)-R)*(cos(p.α)*sin(p.β)p.n[ci]*R^3)
end
function ClusterCondition(theta,p)
    x = theta/2
    1 .-(p.σ*3*18*16*sin(p.α+p.β)/2) .*(sin.(x) ./(2 .*theta)).*((x .*cos.(x) .-sin.(x)) ./(4 .*theta.*theta))
end
function AdlerFactor(state,p)
    (cos(p.β)*(p.n[1]-0.5)-(p.mean[1]-p.mean[2])/p.σ/ContOrder(state,p,1)^2/ContOrder(state,p,2)^2)/sqrt(cos(p.β)^2*(p.n[1]-0.5)^2+sin(p.β)^2/4)
end
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
function ContOrder(state,p::NamedTuple,ci::Integer)
    ContOrder(state[ci],p.n[ci])
end
function ContOrder(theta::Number,n::Number)
    if theta == 0
        1
    else
        4/(theta*n)*sin(theta*n/4)
    end
end
function ContDerivs(du,u,p::NamedTuple,t)
    R=zeros(p.N)
    du[1:end].=0
    for i in 1:p.N
        R[i]=ContOrder(u,p,i)
    end
    for i in 1:p.N
        for j in 1:p.N
            du[i]+=R[j]*p.n[j]*cos(p.α+u[p.N+i]-u[p.N+j])*u[p.N*(i+1)+j]
            du[p.N*(i+1)+j]=-p.ϵ*(u[p.N*(i+1)+j]+R[i]*R[j]*sin(u[p.N+i]-u[p.N+j]+p.β))
            du[p.N+i]-=p.σ*p.n[j]*R[j]*R[i]*u[p.N*(i+1)+j]*sin(u[p.N+i]-u[p.N+j]+p.α)
        end
        du[p.N+i]+=p.mean[i]
        du[i]*=p.σ/p.Var[i]/u[i]*(cos(p.n[i]*u[i]/4)-R[i])
        du[i]+=1
    end
    du
end
function ContOrder(state,params::NamedTuple)
    R = zeros(params.N)
    for i in 1:params.N
        R[i]=ContOrder(state,params,i)
    end
    ContOrder(state,R,params)
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
    sqrt(Order)
end
function BifurcationOrderTuple(sol,params,t::Tuple{Any,Any})
    output1=0
    output2=zeros(params.N)
    for i in 1:length(output2)
        output2[i]=quadgk(x->ContOrder(2 .*sol(x),params,i),t[2]-100,t[2])[1]/(100)
    end
    for i in 1:params.N
        for j in i+1:params.N
            if abs(sol(t[2])[(i+1)*params.N+j])>=0.5
                output1 += 2*params.n[i]*output2[i]*params.n[j]*output2[j]*cos(2*(sol(t[2])[i+params.N]-sol(t[2])[j+params.N]))
            end
        end
        output1 += params.n[i]^2 *output2[i]^2
    end
    if output1<0
        error(sol(t[2]),"\n",output1,"\n",output2)
    end
    sqrt(output1),output2...
end
function BifurcationOrderTupleProperlyTimed(sol,params,t::Tuple{Any,Any};t_anal=1500)
    output1=0
    output2=zeros(params.N)
    for i in 1:length(output2)
        output2[i]=quadgk(x->ContOrder(2 .*sol(x),params,i),t[2]-t_anal,t[2])[1]/(t_anal)
    end
    output1 = quadgk(x->ContOrder(2 .*sol(x),output2,params),t[2]-t_anal,t[2])[1]/(t_anal)
    output1,output2...
end
function TimeaveragedFirstOrder(sol,p,t::Tuple{<:Real,<:Real})
    quadgk(x->ContOrder(sol(x),p),t[1],t[2])[1]/(t[2]-t[1])
end
function ContScanthroughSigma(Initialstate, n1, Measure;Sigmas=collect(0:0.01:10),
t=(0.,500.),α=0.0*pi,Betamult=-53,β=Betamult/100*pi,filename="Ordermatrix",wigglesize=0.01,kwargs...)

    params = CreateContParams(n1...,α=α ,β=β ,σ=Sigmas[1])
    prob=ODEProblem(ContDerivs,Initialstate,(0,0.),params)
    sol=DifferentialEquations.solve(prob)
    Order = Vector{typeof(Measure(sol,params,(0.,0)))}(undef,length(Sigmas))
    state = Initialstate
    for i in 1:length(Sigmas)
        prob = ODEProblem(ContDerivs, state, t, params)
        sol = DifferentialEquations.solve(prob; kwargs...)
        Order[i]=Measure(sol,params,t)
        if i < length(Sigmas)
            params = CreateContParams(n1...,α=α,β=β,σ=Sigmas[i+1])
        end
        wiggle=rand(length(state))
        state#=[1:end#=params.N=#]=# = sol(t[2])#=[1:params.N]=#.+wigglesize.*wiggle#[1:params.N] .-wigglesize/2
        #state[params.N+1:end] = sol(t[2])[params.N+1:end].+wigglesize.*wiggle[params.N+1:end] .-wigglesize/2
        for i in 1:params.N
            if state[i]>10.
                state[i]=sol(0)[i]
                #state[params.N*2+(i-1)*params.N+i]=sol(0)[params.N*2+(i-1)*params.N+i]
            end
        end
    end
    return Order, state
end
function CreateBifDiagramSchar(Initialstate,n1,filepath;sigmas = collect(0:0.01:10),t=(0.,500.),wigglesize=0.01,measure = BifurcationOrderTupleProperlyTimed,kwargs...)
    Orders = zeros((2*length(n1),length(sigmas)+1))
    sigmasreverse = sigmas[end:-1:1]
    Orders[1:end,1]=[n1...,n1...]
    for i in 1:length(n1)
        n=n1[i]
        Order,Finalstate = ContScanthroughSigma(Initialstate, n, measure; β = -53*pi/100,Sigmas=sigmas,t=t,wigglesize=wigglesize,kwargs...)
        Orders[i,2:end]=[e[1] for e in Order]
        Order,Finalstate = ContScanthroughSigma(Finalstate, n, measure; β = -53*pi/100,Sigmas=sigmasreverse,t=t,wigglesize=wigglesize,kwargs...)
        Orders[i+length(n1),2:end]=[e[1] for e in Order[end:-1:1]]
    end
    writedlm(filepath,Orders)
end
function CreateThreeClusterBifDiagramSchar(Initialstate,n1,filepath;sigmas = collect(0:0.01:10),t=(0.,500.),wigglesize=0.01,kwargs...)
    Orders = zeros((2*length(n1),length(sigmas)+2))
    sigmasreverse = sigmas[end:-1:1]
    Orders[1:end,1]=[[e[1] for e in n1]...,[e[1] for e in n1]...]
    Orders[1:end,2]=[[e[2] for e in n1]...,[e[2] for e in n1]...]
    for i in 1:length(n1)
        n=n1[i]
        print(n," ")
        Order,Finalstate = ContScanthroughSigma(Initialstate, n, BifurcationOrderTuple; β = -53*pi/100,Sigmas=sigmas,t=t,wigglesize=wigglesize,kwargs...)
        Orders[i,3:end]=[e[1] for e in Order]
        Order,Finalstate = ContScanthroughSigma(Finalstate, n, BifurcationOrderTuple; β = -53*pi/100,Sigmas=sigmasreverse,t=t,wigglesize=wigglesize,kwargs...)
        Orders[i+length(n1),3:end]=[e[1] for e in Order[end:-1:1]]
        #=
        for i in 1:1
            plot!(sigmas,[e[i] for e in Order],c=get(ColorSchemes.bamako,n,(0.5,1.)),style=:dash,label="")
        end
        =#
    end
    writedlm(filepath,Orders)
end
function PlotBifDiagram(Orders,sigmas;skips=1,downsweep=true,cutoff=0.01,kwargs...)
    n=size(Orders)[1]
    if downsweep
        n = Int(size(Orders)[1]/2)
    end
    for i in 1:n
        previndex = 1
        for j in (skips+1):length(sigmas)
            if abs(Orders[i,j+skips]-Orders[i,j+skips-1])>cutoff
                plot!(sigmas[previndex:j-1],Orders[i,previndex+skips:j+skips-1],c=get(ColorSchemes.bamako,Orders[i,1],(0.5,1.0)),label="",kwargs...)
                previndex = j
            end
        end
        plot!(sigmas[previndex:end],Orders[i,previndex+skips:end],c=get(ColorSchemes.bamako,Orders[i,1],(0.5,1.0)),label=string("n_1 = ",Orders[i,1]),kwargs...)
        previndex = 1
        if downsweep
            plot!(sigmas[1:end],Orders[i+n,skips+1:end],c=get(ColorSchemes.bamako,Orders[i+n,1],(0.5,1.0)),label="",style=:dash,kwargs...)
        end
    end
    plot!(xlabel = "Sigma", ylabel="2nd Orderparameter")
end
seed = 31#138513#parse(Int64, ARGS[1])
Random.seed!(seed)
initstate = rand(8)
initstate[3:4].*=2*pi
println(initstate)

function CreateBifDiagramScharUpdated(Initialstate,n1,filepath,sigmastart;sigmas = collect(0:0.01:10),t=(0.,500.),wigglesize=0.01,measure = BifurcationOrderTupleProperlyTimed,kwargs...)
    Orders = zeros((2*length(n1),length(sigmas)+1))
    sigma1 = collect(sigmastart:sigmas[2]-sigmas[1]:sigmas[end])
    sigma2 = collect(sigmas[1]:sigmas[2]-sigmas[1]:sigmastart)
    sigmasreverse = sigmas[end:-1:1]
    Orders[1:end,1]=[n1...,n1...]
    for i in 1:length(n1)
        n=n1[i]
        Order,Finalstate = ContScanthroughSigma(Initialstate, n, measure; β = -53*pi/100,Sigmas=sigma1,t=t,wigglesize=wigglesize,kwargs...)
        Orders[i,1+findfirst(isequal(sigma1[1]),sigmas):end]=[e[1] for e in Order]
        Order,Finalstate = ContScanthroughSigma(Finalstate, n, measure; β = -53*pi/100,Sigmas=sigmasreverse,t=t,wigglesize=wigglesize,kwargs...)
        Orders[i+length(n1),2:end]=[e[1] for e in Order[end:-1:1]]
        Order,Finalstate = ContScanthroughSigma(Initialstate, n, measure; β = -53*pi/100,Sigmas=sigma2[end:-1:1],t=t,wigglesize=wigglesize,kwargs...)
        Orders[i,2:1+findfirst(isequal(sigma1[1]),sigmas)]=[e[1] for e in Order]

    end
    writedlm(filepath,Orders)
end
CreateBifDiagramScharUpdated(initstate,collect(0.5:0.05:0.95),"./Data/Singleruns/Orders.txt",0.5,sigmas=collect(0.:0.01:5),t=(0.,2000.),dtmax=10.10,wigglesize=0.01,measure=BifurcationOrderTupleProperlyTimed)
#CreateBifDiagramSchar(initstate,collect(0.5:0.05:0.95),"./Data/Singleruns/OrdersShort.txt",t=(0.,1500.),sigmas=collect(0.6:-0.01:0),dtmax=10,wigglesize=0.01,measure=BifurcationOrderTupleProperlyTimed)
