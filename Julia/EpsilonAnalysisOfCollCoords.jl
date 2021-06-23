using Random
using Roots
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
    params = (omegas = omegas, indices = indices, Σ = S, mean = m, α=α, β=β, σ=σ, N=n, ϵ=ϵ)
end
function WeightedSumOverModCosPairs(minuend,subtrahend,weights)
    sinm = sin.(minuend)
    cosm = cos.(minuend)
    sins = sin.(subtrahend)
    coss = cos.(subtrahend)
    WeightedSumOverModCosPairs(sinm,cosm,sins,coss,weights)
end
function WeightedSumOverModCosPairs(sinm,cosm,sins,coss,weight)
    output=0
    intermediate = 0
    for i in 1:length(sinm)
        for j in 1:length(sins)
            intermediate+=(cosm[i]coss[j]+sinm[i]sins[j])*(weight[i]-weight[j])
        end
        output+=weight[i]*intermediate
        intermediate=0
    end
    output
end
function WeightedSumOverSinePairs(minuend,subtrahend,weights)
    sinm = sin.(minuend)
    cosm = cos.(minuend)
    sins = sin.(subtrahend)
    coss = cos.(subtrahend)
    WeightedSumOverSinePairs(sinm,cosm,sins,coss,weights)
end
function WeightedSumOverSinePairs(sinm,cosm,sins,coss,weight)
    output=0
    intermediate = 0
    for i in 1:length(sinm)
        for j in 1:length(sins)
            intermediate+=sinm[i]coss[j]-cosm[i]sins[j]
        end
        output+=weight[i]*intermediate
        intermediate=0
    end
    output
end
function SumOverSinePairs(minuend,subtrahend)
    sinm = sin.(minuend)
    cosm = cos.(minuend)
    sins = sin.(subtrahend)
    coss = cos.(subtrahend)
    SumOverSinePairs(sinm,cosm,sins,coss)
end
function SumOverSinePairs(sinm,cosm,sins,coss)
    output=0
    for i in 1:length(sinm)
        for j in 1:length(sins)
            output+=sinm[i]coss[j]-cosm[i]sins[j]
        end
    end
    output
end
function Kappa0(theta0,params,ci;
    sin0 = sin.(theta0.*params.omegas[params.indices[ci]]),
    cos0 = cos.(theta0.*params.omegas[params.indices[ci]]),
    sinb = sin.(theta0.*params.omegas[params.indices[ci]].+params.β),
    cosb = cos.(theta0.*params.omegas[params.indices[ci]].+params.β))
    output = SumOverSinePairs(sinb,cosb,sin0,cos0)
    N = length(params.indices[ci])
    -output/(N*N)
end
function Theta0(theta0,params,ci;
    sin0 = sin.(theta0.*params.omegas[params.indices[ci]]),
    cos0 = cos.(theta0.*params.omegas[params.indices[ci]]),
    sina = sin.(theta0.*params.omegas[params.indices[ci]].+params.α),
    cosa = cos.(theta0.*params.omegas[params.indices[ci]].+params.α))
    kappa0 = Kappa0(theta0,params,ci,sin0=sin0,cos0=cos0)
    output = WeightedSumOverSinePairs(sina,cosa,sin0,cos0,params.omegas[params.indices[ci]].-params.mean[ci])
    output*=kappa0/length(params.indices[ci])/length(omegas)*params.σ/params.Σ[ci]
    output-=1
end
function Omega0(theta0,params,ci;
    sin0 = sin.(theta0.*params.omegas[params.indices[ci]]),
    cos0 = cos.(theta0.*params.omegas[params.indices[ci]]),
    sina = sin.(theta0.*params.omegas[params.indices[ci]].+params.α),
    cosa = cos.(theta0.*params.omegas[params.indices[ci]].+params.α))
    output = SumOverSinePairs(sina,cosa,sin0,cos0)
    output*=Kappa0(theta0,params,ci,sin0=sin0,cos0=cos0)/length(params.indices[ci])
    params.mean[ci]-params.σ/length(params.omegas)*output
end
function Omegaprime(theta0s,params)
    output1=Omega0(theta0s[1],params,1)
    output2=Omega0(theta0s[2],params,2)
    output1-output2
end
function FindTheta0(params,ci)
    find_zero(x -> Theta0(x,params,ci),0)
end
Random.seed!(0)
omegas = sort(0.5*rand(50)-0.25.*ones(50))
p = CreateParams(omegas,[collect(1:27),collect(28:50)],α=0.0*pi,β=-0.53*pi,σ=3.5)
T0 = FindTheta0(p,1)
thetas=collect(-50:0.01:50)
rhs=[Theta0(x,p,1) for x in thetas]
