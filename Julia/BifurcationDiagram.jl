using DifferentialEquations
using Random
using Statistics
using DelimitedFiles
using Plots
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
function ContDerivsOLD(du,u,p::NamedTuple,t)
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
function CreateBifDiagram(Initialstate,n1,filepath)
    sigmas = collect(0.:0.01:10)
    colors = ["blue","green","red"]
    labels=["Full System","First Cluster","Second Cluster"]
    Order,Finalstate = ContScanthroughSigma(Initialstate, n1, BifurcationOrderTuple, β = -53*pi/100,Sigmas=sigmas)
    for i in 1:length(Order[1])
        plot!(sigmas,[e[i] for e in Order],c=colors[i],style=:solid,label=labels[i])
    end
    sigmas=collect(10:-0.01:0)
    Order,Finalstate = ContScanthroughSigma(Finalstate, n1, BifurcationOrderTuple, β = -53*pi/100,Sigmas=sigmas)
    for i in 1:length(Order[1])
        plot!(sigmas,[e[i] for e in Order],c=colors[i],style=:dash,label="")
    end
    plot!(xlabel="Sigma",ylabel="Orderparameter")
    savefig(filepath)
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
function PlotBifDiagram(Orders,sigmas;skips=1)
    n = Int(size(Orders)[1]/2)
    for i in 1:n
        previndex = 1
        for j in (skips+1):length(sigmas)
            if abs(Orders[i,j+skips]-Orders[i,j+skips-1])>0.01
                plot!(sigmas[previndex:j-1],Orders[i,previndex+skips:j+skips-1],c=get(ColorSchemes.bamako,Orders[i,1],(0.5,1.0)),label="")
                previndex = j
            end
        end
        plot!(sigmas[previndex:end],Orders[i,previndex+skips:end],c=get(ColorSchemes.bamako,Orders[i,1],(0.5,1.0)),label=string("n_1 = ",Orders[i,1]))
        previndex = 1
        plot!(sigmas[1:end],Orders[i+n,skips+1:end],c=get(ColorSchemes.bamako,Orders[i+n,1],(0.5,1.0)),label="",style=:dash)
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
