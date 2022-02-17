using LinearAlgebra
using Plots
using Random
using Roots
using Serialization

function ThetaCondition(theta,omegas,sigma,a=0,b=-0.53*pi)
    N = length(omegas)
    output = sum(omegas .* omegas)/N
    cw = cos.(2*theta .* omegas)
    sw = sin.(2*theta .* omegas)
    ca = cos(a+b)
    sa = sin(a+b)
    output2 = 0
    for i in 1:N
        runner = 0
        for j in 1:N
            runner += cw[i]*(ca*cw[j]+sa*sw[j])-sw[i]*(sa*cw[j]-ca*sw[j])
        end
        output2 += omegas[i]*runner
    end
    output-output2*sigma/N/N/2
end
function MatrixA(phases,sigma,α=0,β=-0.53*pi)
    N = length(phases)
    A = zeros((N,N))
    for i in 1:N
        for j in 1:N
            if i!=j
                A[i,j]=-sigma/N*sin(phases[i]-phases[j]+β)*cos(phases[i]-phases[j]+α)
            end
        end
        A[i,i]=-sum(A[i,1:end])
    end
    A
end
function MatrixB(phases,n,sigma,α=0,β=-0.53*pi)
    N = length(phases)
    B = zeros((N,N))
    for j in 1:N
        B[n,j]=-sigma/N*sin(phases[n]-phases[j]+α)
    end
    B
end
function MatrixC(phases,n,epsilon,α=0,β=-0.53*pi)
    N = length(phases)
    C = zeros((N,N))
    for i in 1:N
        if i!=n
            C[i,n] = -epsilon*cos(phases[n]-phases[i]+β)
            C[i,i] = epsilon*cos(phases[n]-phases[i]+β)
        end
    end
    C
end
function Jacobian(phases,sigma,epsilon,α=0,β=-0.53*pi)
    N = length(phases)
    A = MatrixA(phases,sigma,α,β)
    B = MatrixB(phases,1,sigma,α,β)
    C = MatrixC(phases,1,epsilon,α,β)
    for i in 2:N
        B = [B MatrixB(phases,i,sigma,α,β)]
        C = [C; MatrixC(phases,i,epsilon,α,β)]
    end
    Jacobian = [A B; C -epsilon*I]
end

Random.seed!(10)
omegas = sort!(0.5*rand(50) .- 0.25)
n1 = collect(1:50)
sigmas = collect(0:0.5:5)
stability = zeros(Float64,length(n1),length(sigmas))
α=0
β=-0.53*pi
for n1 in n1[1:end-1]
    omegas1 = omegas[1:n1]
    omegas2 = omegas[n1+1:end]
    println(length(omegas1),length(omegas2))
    for i in 1:length(sigmas)
        sigma = sigmas[i]
        theta1=0
        theta2=0
        try
            theta1 = find_zero(x -> ThetaCondition(x,omegas1,sigma,α,β),0)
            if theta1<-0 || theta1>10
                theta1 = missing
            end
        catch
            theta1 = missing
        end
        try
            theta2 = find_zero(x -> ThetaCondition(x,omegas2,sigma,α,β),0)
            if theta2<-0 || theta2>10
                theta2 = missing
            end
        catch
            theta2 = missing
        end
        if ismissing(theta1) || ismissing(theta2)
            stability[n1,i] = NaN
        else
            result1 = sort(real.(eigvals(Jacobian(theta1 .* omegas1,sigma*n1/50,0.01,α,β))))
            result2 = sort(real.(eigvals(Jacobian(theta2 .* omegas2,sigma*(1-n1/50),0.01,α,β))))
            r1 = result1[end]<10^-10 ? result1[end-1] : result1[end]
            r2 = result2[end]<10^-10 ? result2[end-1] : result2[end]
            stability[n1,i] = r1<r2 ? r2 : r1
        end
    end
end
for i in 1:length(sigmas)
    sigma=sigmas[i]
    theta=0
    try
        theta = find_zero(x -> ThetaCondition(x,omegas,sigma,α,β),0)
        if theta<0 || theta>10
            theta = missing
        end
    catch
        theta = missing
    end
    if ismissing(theta)
        stability[50,i]=NaN
    else
        result2 = sort(real.(eigvals(Jacobian(theta .* omegas,sigma,0.01,α,β))))
        stability[50,i]=result2[end]<10^-10 ? result2[end-1] : result2[end]
    end
end
heatmap(sigmas,n1,stability,xlabel="Sigma",ylabel="Oscillators in Cluster 1")
open("./Stability_negativeThetas.serialjl","w") do io
    serialize(io,stability)
end
#=
alphas = range(0,0.5*pi,length=25)
betas = range(0,2*pi,length=100)
for i in eachindex(alphas)
    for j in eachindex(betas)
       stability[i,j] = sort(real.(eigvals(Jacobian(phases,1,0.01,alphas[i],betas[j]))))[end]<0.000001 ? 1 : 0
   end
end
heatmap(stability)
=#
#=
Random.seed!(0)
omegas = sort!(0.5 .* rand(50) .- 0.25)
thetas = collect(0:0.02:2)
sigmas = collect(0:0.01:1)
stabilities = zeros(length(sigmas),length(thetas))
for i in 1:length(sigmas)
    for j in 1:length(thetas)
        results = sort(real.(eigvals(Jacobian(thetas[j] .* omegas,sigmas[i],0.01,0.3*pi,-0.53*pi))))
        stabilities[i,j] = results[end]
    end
end
=#
#=
LLE=zeros(length(thetas))
for i in eachindex(LLE)
    results = sort(real.(eigvals(Jacobian(thetas[i] .* omegas,0.1,0.01,0.3*pi,-0.53*pi))))
    LLE[i]=results[end]#+1 ≈ 1 ? results[end-1] : results[end]
end
=#
