using DifferentialEquations
using Random
using Statistics
using DelimitedFiles
using Plots
using QuadGK
using Roots

#Similar to BifurcationDiagram, but for the new Ansatz, changes: New Derivative, new Parameters, new OrderTuple.
function SteadyClusterStateCondition(state,p,ci)
    SteadyClusterStateCondition(state[ci],p,ci)
end
function SteadyClusterStateCondition(theta::Number,p,ci)
    R = ContOrder(theta,p.n[ci])
    1-p.σ/theta/p.Var[ci]*(cos(theta*p.n[ci]/4)-R)*(cos(p.α)*sin(p.β)p.n[ci]*R^3)
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
    OPrime = min(OmegaPrime(SteadyStateOrders(params),params)...)
    params = (N=N, Var = var, mean = mean, OPrime = OPrime, α=α, β=β, σ=σ, ϵ=ϵ,n=[n...])
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
function SteadyStateTheta(theta::Number, p, ci)
    R2 = ContOrder(2*theta,p.n[ci])
    1-p.σ/2/p.Var[ci]/theta*(cos(theta*p.n[ci]/2)-R2)*p.n[ci]*R2*sin(p.β)
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
        R[i]=ContOrder(2 .*u,p,i)
    end
    du[1] = 1-p.σ/2/p.Var[1]/u[1]*(cos(u[1]*p.n[1]/2)-R[1])*(p.n[1]*R[1]*sin(p.β)-p.ϵ*p.n[2]/p.OPrime*R[2]*cos(2*p.OPrime*t+2*u[3]+p.β))
    du[2] = 1-p.σ/2/p.Var[2]/u[2]*(cos(u[2]*p.n[2]/2)-R[2])*(p.n[2]*R[2]*sin(p.β)+p.ϵ*p.n[1]/p.OPrime*R[1]*cos(-2*p.OPrime*t-2*u[3]+p.β))
    du[3] = p.mean[1]-p.mean[2]-p.OPrime+p.σ*cos(p.β)/2*(p.n[1]*(1-R[1]^2)-p.n[2]*(1-R[2]^2))-p.σ*p.ϵ/2/p.OPrime*(R[1]*R[2]*(p.n[2]*sin(2*p.OPrime*t+2*u[3]+p.β)+p.n[1]*sin(p.β-2*p.OPrime*t-2*u[3]))-sin(p.β))
    du
end
function ContOrder(state,params::NamedTuple,t)
    R = zeros(params.N)
    for i in 1:params.N
        R[i]=ContOrder(state,params,i)
    end
    ContOrder(state,R,params,t)
end
function ContOrder(state,orders,params::NamedTuple,t)
    R = orders
    Order = 0
    for i in 1:params.N
        Order+=params.n[i]*params.n[i]*R[i]*R[i]
    end
    Order+=2*params.n[1]*params.n[2]*R[1]*R[2]*cos(state[3]+params.OPrime*t)
    sqrt(Order)
end
function BifurcationOrderTuple(sol,params,t::Tuple{Any,Any})
    output1=0
    output2=zeros(params.N)
    for i in 1:length(output2)
        output2[i]=quadgk(x->ContOrder(2 .*sol(x),params,i),t[1],t[2])[1]/(t[2]-t[1])
    end
    #=
    for i in 1:params.N
        for j in i+1:params.N
            if abs(sol(t[2])[(i+1)*params.N+j])>=0.5
                output1 += 2*params.n[i]*output2[i]*params.n[j]*output2[j]*cos(2*(sol(t[2])[i+params.N]-sol(t[2])[j+params.N]))
            end
        end
        output1 += params.n[i]^2 *output2[i]^2
    end
    =#
    output1 = quadgk(t -> ContOrder(2 .*sol(t),params,t),t[1],t[2])[1]/(t[2]-t[1])
    if output1<0
        error(sol(t[2]),"\n",output1,"\n",output2)
    end
    output1,output2...
end
function TimeaveragedFirstOrder(sol,p,t::Tuple{<:Real,<:Real})
    quadgk(x->ContOrder(sol(x),p),t[1],t[2])[1]/(t[2]-t[1])
end
function ContScanthroughSigma(Initialstate, n1, Measure;Sigmas=collect(0:0.01:10),
t=(0.,500.),α=0.0*pi,Betamult=-53,β=Betamult/100*pi,filename="Ordermatrix",wigglesize=0.01,t_anal=(t[2]-100,t[2]),kwargs...)

    params = CreateContParams(n1...,α=α ,β=β ,σ=Sigmas[1])
    prob=ODEProblem(ContDerivs,Initialstate,(0,0.),params)
    sol=DifferentialEquations.solve(prob)
    Order = Vector{typeof(Measure(sol,params,(0.,0)))}(undef,length(Sigmas))
    state = Initialstate
    for i in 1:length(Sigmas)
        prob = ODEProblem(ContDerivs, state, t, params)
        sol = DifferentialEquations.solve(prob; kwargs...)
        Order[i]=Measure(sol,params,t_anal)
        if i < length(Sigmas)
            params = CreateContParams(n1...,α=α,β=β,σ=Sigmas[i+1])
        end
        #=
        if rem(Sigmas[i],1)==0
            writedlm(string("./",n1*50,"_",seed,"_",Sigmas[i],"_",filename,".txt"),sol)
        end
        =#
        state = sol(100).+wigglesize.*rand(length(state))
        for i in 1:params.N
            if state[i]>10.
                state[i]=sol(0)[i]
            end
        end
    end
    return Order, state
end
function CreateBifDiagramData(Initialstate,n1,filepath;t=(0.,500.))
    sigmas = collect(0.:0.01:10)
    Order,Finalstate = ContScanthroughSigma(Initialstate, n1, BifurcationOrderTuple, β = -53*pi/100,Sigmas=sigmas,t=t)
    writedlm(filepath[1],Order)
    sigmas=collect(10:-0.01:0)
    Order,Finalstate = ContScanthroughSigma(Finalstate, n1, BifurcationOrderTuple, β = -53*pi/100,Sigmas=sigmas,t=t)
    writedlm(filepath[2],Order)
end
function CreateBifDiagramSchar(Initialstate,n1,filepath;sigmas = collect(0:0.01:10),t=(0.,500.),wigglesize=0.01,kwargs...)
    Orders = zeros((2*length(n1),length(sigmas)+1))
    sigmasreverse = sigmas[end:-1:1]
    Orders[1:end,1]=[n1...,n1...]
    for i in 1:length(n1)
        n=n1[i]
        Order,Finalstate = ContScanthroughSigma(Initialstate, n, BifurcationOrderTuple; β = -53*pi/100,Sigmas=sigmas,t=t,wigglesize=wigglesize,kwargs...)
        Orders[i,2:end]=[e[1] for e in Order]
        Order,Finalstate = ContScanthroughSigma(Finalstate, n, BifurcationOrderTuple; β = -53*pi/100,Sigmas=sigmasreverse,t=t,wigglesize=wigglesize,kwargs...)
        Orders[i+length(n1),2:end]=[e[1] for e in Order[end:-1:1]]
        #=
        for i in 1:1
            plot!(sigmas,[e[i] for e in Order],c=get(ColorSchemes.bamako,n,(0.5,1.)),style=:dash,label="")
        end
        =#
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
function CreateBifDiagramZoomed(Initialstate,n1,filepath)
    sigmas = collect(0.:0.01:10)
    labels=["Full System","First Cluster","Second Cluster"]
    Order,Finalstate = ContScanthroughSigma(Initialstate, n1, BifurcationOrderTuple, β = -53*pi/100,Sigmas=sigmas)
    for i in 1:length(Order[1])
        plot!(sigmas,[e[i] for e in Order],style=:solid,label=labels[i])
    end
    sigmas=collect(10:-0.01:0)
    Order,Finalstate = ContScanthroughSigma(Finalstate, n1, BifurcationOrderTuple, β = -53*pi/100,Sigmas=sigmas)
    for i in 1:length(Order[1])
        plot!(sigmas,[e[i] for e in Order],style=:dash,label=labels[i])
    end
    plot!(xlabel="Sigma",ylabel="2nd Orderparameter",ylim = [0.99,1])
    savefig(filepath)
end
function PlotBifDiagram(Orders,sigmas;skips=1,downsweep=true,cutoff=0.01)
    n=size(Orders)[1]
    if downsweep
        n = Int(size(Orders)[1]/2)
    end
    for i in 1:n
        previndex = 1
        for j in (skips+1):length(sigmas)
            if abs(Orders[i,j+skips]-Orders[i,j+skips-1])>cutoff
                plot!(sigmas[previndex:j-1],Orders[i,previndex+skips:j+skips-1],c=get(ColorSchemes.bamako,Orders[i,1],(0.5,1.0)),label="")
                previndex = j
            end
        end
        plot!(sigmas[previndex:end],Orders[i,previndex+skips:end],c=get(ColorSchemes.bamako,Orders[i,1],(0.5,1.0)),label=string("n_1 = ",Orders[i,1]))
        previndex = 1
        if downsweep
            plot!(sigmas[1:end],Orders[i+n,skips+1:end],c=get(ColorSchemes.bamako,Orders[i+n,1],(0.5,1.0)),label="",style=:dash)
        end
    end
    plot!(xlabel = "Sigma", ylabel="2nd Orderparameter")
end
seed = 138513#parse(Int64, ARGS[1])
n1mult = 27#parse(Int64,ARGS[2])
Random.seed!(seed)
initstate = rand(8)
#n1=n1mult/50
#Order,Finalstate=ContScanthroughSigma(initstate,n1,TimeaveragedFirstOrder,Sigmas=collect(0:0.01:1),filename="test2")
#=
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
=#
