function createsins(u, params; func = sin, extra = 0)
    sins = zeros(length(params.omegas))
    for i in 1:params.N
        for k in params.indices[i]
            sins[k] = func(u[i]*(params.omegas[k]-params.mean[i])+u[params.N+i]+extra)
        end
    end
    return sins
end

function derivativestheta(du, u,  params, t; sina=createsins(u,params, extra = params.α), cosa = createsins(u,params, func = cos, extra = params.α), sin0 = createsins(u,params), cos0 = createsins(u,params, func = cos))
    du[1:params.N].=0.0
    for i in 1:params.N#Index des clusters von theta
        for j in 1:params.N#Index der anderen Cluster
            for k in params.indices[i]#Oszillatorindex im Thetacluster
                for l in params.indices[j]#Oszillatorindex im anderen Cluster
                    du[i]-=u[2*params.N+params.N*(i-1)+j]*(params.omegas[k]-params.mean[i])*(sina[k]*cos0[l]-cosa[k]*sin0[l])
                end
            end
        end
        du[i]*=params.σ/length(params.omegas)/params.Σ[i]/size(params.indices[i])[1]
    end
    du[1:params.N].+=1.0
end

function derivativesf(du, u,  params, t; sina=createsins(u,params, extra = params.α), cosa = createsins(u,params, func = cos, extra = params.α), sin0 = createsins(u,params), cos0 = createsins(u,params, func = cos))
    du[params.N+1:2*params.N].=0.0
    for i in 1:params.N#Index des clusters von F
        for j in 1:params.N#Index der anderen Cluster
            for k in params.indices[i]#Oszillatorindex im F-Cluster
                for l in params.indices[j]#Oszillatorindex im anderen Cluster
                    du[params.N+i]-=u[2*params.N+params.N*(i-1)+j]*(sina[k]*cos0[l]-cosa[k]*sin0[l])
                end
            end
        end
        du[params.N+i]*=params.σ/length(params.omegas)/size(params.indices[i])[1]
    end
    du[params.N+1:2*params.N].+=params.mean
end

function derivativeskappa(du, u,  params, t; sinb=createsins(u,params, extra = params.β), cosb = createsins(u,params, func = cos, extra = params.β), sin0 = createsins(u,params), cos0 = createsins(u,params, func = cos))
    du[2*params.N+1:end].=0.0
    for i in 1:params.N#Index des ersten clusters
        for j in 1:params.N#Index des zweiten Cluster
            for k in params.indices[i]#Oszillatorindex im ersten Clsuter
                for l in params.indices[j]#Oszillatorindex im anderen Cluster
                    du[2*params.N+params.N*(i-1)+j]+=(sinb[k]*cos0[l]-cosb[k]*sin0[l])
                end
            end
            du[2*params.N+params.N*(i-1)+j]/=(size(params.indices[i])[1]*size(params.indices[j])[1])
        end
    end
    du[2*params.N+1:end].+=u[2*params.N+1:end]
    du[2*params.N+1:end].*=-params.ϵ
end

function derivativesCollCoords(du, u,  params, t)
    sina = createsins(u,params, extra = params.α)
    cosa = createsins(u,params, func = cos, extra = params.α)
    sin0 = createsins(u,params)
    cos0 = createsins(u,params, func = cos)
    derivativestheta(du, u, params, t, sina=sina, sin0=sin0, cosa=cosa, cos0=cos0)
    derivativesf(du, u, params, t, sina=sina, sin0=sin0, cosa=cosa, cos0=cos0)
    derivativeskappa(du, u, params, t, sin0=sin0, cos0=cos0)
end

function ThetastoPhis(u,params)
    Phis = zeros(size(params.omegas))
    for i in eachindex(params.indices)
        for j in params.indices[i]
            Phis[j]=u[i]*(params.omegas[j]-params.mean[i])+u[params.N+i]
        end
    end
    return Phis
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
    params = (omegas = omegas, indices = indices, Σ = S, mean = m, α=α, β=β, σ=σ, N=n, ϵ=ϵ)
end

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
