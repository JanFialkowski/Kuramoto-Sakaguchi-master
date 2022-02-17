using Random, DifferentialEquations, Statistics, DelimitedFiles, Serialization

struct Kuramotosystem
    state::Array{Float64,1}
    N::Int64
    α::Float64
    β::Float64
    ϵ::Float64
    σ::Float64
    omegas::Array{Float64,1}
end

function ConstructStatefromPhis(phis;β=pi/2)
    N=length(phis)
    u0 = [phis;ones(N*N)]
    for i in N+1:length(u0)
        u0[i] = -sin(u0[div(i-N-1,N)+1]-u0[rem(i-N-1,N)+1]+β)
    end
    return u0
end

function Kuramotosystem(phis; N=length(phis), omegas = [0. for i in 1:N], α=0., β = pi/2, σ = 2.5, ϵ = 0.01)
    state = ConstructStatefromPhis(phis,β=β)
    return Kuramotosystem(state,N,α,β,ϵ,σ,omegas)
end

function Derivativesphi(dstate, system)
    dstate[1:system.N] .= 0
    sina = sin.(system.state[1:system.N] .+system.α)
    cosa = cos.(system.state[1:system.N] .+system.α)
    sin0 = sin.(system.state[1:system.N])
    cos0 = cos.(system.state[1:system.N])
    for i = 1:system.N
        for j = 1:system.N
            dstate[i] -= system.state[system.N*i+j] * (sina[i]*cos0[j]-cosa[i]*sin0[j])
        end
        dstate[i] *= system.σ / system.N
        dstate[i] += system.omegas[i]
    end
end

function Derivativeskappa(dstate, system)
    dstate[system.N+1:end] .= 0
    sinb = sin.(system.state[1:system.N] .+system.β)
    cosb = cos.(system.state[1:system.N] .+system.β)
    sin0 = sin.(system.state[1:system.N])
    cos0 = cos.(system.state[1:system.N])
    for i in 1:system.N
        for j in 1:system.N
            dstate[system.N*i+j] = -system.ϵ*(system.state[system.N*i+j]+(sinb[i]*cos0[j]-cosb[i]*sin0[j]))
        end
    end
end

function order(system::Kuramotosystem)
    Phis = system.state[1:system.N]
    order(Phis)
end

function order(Phis::Array{Float64,1})
    r = 0.0 + 0im
    for i = 1:length(Phis)
        r += exp(1im * Phis[i])
    end
    r /= length(Phis)
    r = abs(r)
end

function order(Phis::Array{Float64,2})
    r = zeros(typeof(1.0 + 1im), length(Phis[1, :]))
    for i = 1:length(Phis[1, :])
        r[i] = order(Phis[:,i])
    end
    r
end

function CustomRK4(system, tspan; dt = 0.1, save_indices=1:system.N, save_times=collect(tspan[1]:dt:tspan[2]))
    sort!(save_times)
    integrated_times = collect(tspan[1]:dt:tspan[2])
    for i in save_times
        if i∉integrated_times
            error("One or more of the specified saved timesteps is not actually calculated.")
        end
    end
    dstate = zeros(length(system.state))
    phistep = zeros(system.N)
    phibetween = zeros(system.N)
    save_times_iterator = 2
    if tspan[1]∉save_times
        pushfirst!(save_times,tspan[1])
    end
    saved_stuff = (t = save_times,
    u = zeros(length(system.state[save_indices]), length(save_times)),)
    saved_stuff.u[:, 1] .= system.state[save_indices]
    for i = 2:length(integrated_times)
        Derivativesphi(dstate, system)
        phistep = dstate[1:system.N] #phistep = k1
        phibetween = dt * phistep / 2
        system.state[1:system.N] += phibetween #u+k1/2
        Derivativesphi(dstate, system)
        phistep += 2 * dstate[1:system.N] #phistep = k1+2k2
        system.state[1:system.N] -= phibetween #u
        phibetween = dt * dstate[1:system.N] / 2 #k2/2*dt
        system.state[1:system.N] += phibetween #u+k2/2*dt
        Derivativesphi(dstate, system)
        phistep += 2 * dstate[1:system.N] #phistep = k1+2k2+2k3
        system.state[1:system.N] -= phibetween #u
        phibetween = dt * dstate[1:system.N] #k3*dt
        system.state[1:system.N] += phibetween #u+k3*dt
        Derivativesphi(dstate, system)
        phistep += dstate[1:system.N] #phistep = k1+2k2+2k3+k4
        system.state[1:system.N] -= phibetween #u0
        Derivativeskappa(dstate, system)
        system.state[1:system.N] += dt * phistep / 6
        system.state[system.N+1:end] += dt * dstate[system.N+1:end]
        if integrated_times[i]==save_times[save_times_iterator]
            saved_stuff.u[:, save_times_iterator] .= system.state[save_indices]
            if save_times_iterator < length(save_times)
                save_times_iterator+=1
            end
        end
        #if rem(100 * (i-1), length(tspan[1]:dt:tspan[2])-1) == 0
        #    println(100 * (i-1) / (length(tspan[1]:dt:tspan[2])-1))
        #end
    end
    return saved_stuff
end

function SetupState(seed,Oszis,n)
    Random.seed!(seed)
    omegas = sort!(0.5*rand(Oszis)-0.25 .*ones(Oszis))
    state= zeros(Oszis+Oszis*Oszis)
    for i=1:Oszis, j=1:Oszis
        if i<=n
            state[i]=0.3*pi
            if j<=n
                state[Oszis*i+j]=-sin(-0.53*pi)
            end
        elseif i>n
            if j>n
                state[Oszis*i+j]=-sin(-0.53*pi)
            end
        end

    end
    return state, omegas
end

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

function CreateParams2(omegas, indices; α=0.0*pi, β=-0.53*pi, σ=3.5, ϵ=0.01)
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

function CalcOrderMatrix(Ordermatrix, sol, Oszis; t1 = 3000, t2 = 10000)
    for k = 1:Oszis, l = 1:Oszis
        if abs(
            (sol(t2)[k] - sol(t1)[k]) / (t2 - t1) -
            (sol(t2)[l] - sol(t1)[l]) / (t2 - t1),
        ) <= 0.001
            Ordermatrix[k, l] = 1.0
        else
            Ordermatrix[k, l] = 0
        end
    end
end
function Analyze(sol,n,Orderovertime,Orders,phis,i)
    if n != 1000
        Orderovertime = order(2 .* sol.u[1:n,:])
        Orders[2,i]=mean(abs.(Orderovertime[2:end]))
        Orderovertime = order(2 .* sol.u[n+1:end,:])
        Orders[3,i]=mean(abs.(Orderovertime[2:end]))
        Orderovertime = order(2 .* sol.u)
        Orders[1,i]=mean(abs.(Orderovertime[2:end]))
        phis[1,i]=sol.u[1,2]
        phis[2,i]=sol.u[end,2]
        phis[3,i]=sol.u[1,end]
        phis[4,i]=sol.u[end,end]
        Orderovertime
    else
        Orderovertime = order(2 .* sol.u[1:n,:])
        Orders[2,i]=mean(abs.(Orderovertime[2:end]))
        Orders[3,i]=mean(abs.(Orderovertime[2:end]))
        Orders[1,i]=mean(abs.(Orderovertime[2:end]))
        phis[1,i]=sol.u[1,2]
        phis[2,i]=sol.u[end,2]
        phis[3,i]=sol.u[1,end]
        phis[4,i]=sol.u[end,end]
        Orderovertime
    end
end
function Run(System,tspan;Analspan=1000.)
    sol = CustomRK4(System,tspan,dt=0.01,save_times = collect(tspan[2]-Analspan:0.1:tspan[2]))
end
function SetupSystem(n,sigma)
    state,omegas = SetupState(0,1000,n)
    Kuramotosystem(state,1000,0,-0.53*pi,0.01,sigma,omegas)
end
function Dothething(
    n,
    sigmas = collect(1.0:0.1:5.5);
    state = nothing,
    Extrastring="",
    save_end=true,
    save_phis=false,
    Analspan=1000.,
    tspan = (0.,1100.)
    )

    Orderovertime = zeros(Int(Analspan)*10+2)
    Orders = zeros((3,length(sigmas)))
    phis = zeros((4,length(sigmas)))
    System = SetupSystem(n,sigmas[1])
        if state != nothing
            System.state[1:end] .= state
        end
    for i in 1:length(sigmas)
        System = Kuramotosystem(System.state,1000,0,-0.53*pi,0.01,sigmas[i],System.omegas)
        sol = Run(System,tspan;Analspan=Analspan)
        Orderovertime=Analyze(sol,n,Orderovertime,Orders,phis,i)
        if i>2 && abs(Orders[1,i]-Orders[1,i-1])>2*abs(Orders[1,i-1]-Orders[1,i-2])
            sol = Run(System,(tpsan[1],10*tspan[2]);Analspan=Analspan)
            Orderovertime=Analyze(sol,n,Orderovertime,Orders,phis,i)
        end
        writedlm(string("./",Extrastring,"Timeline_",sigmas[i],".txt"),Orderovertime)
        print(sigmas[i])
        if save_phis
            writedlm(string("./",Extrastring,"PhiTimeline_",sigmas[i],".txt"),sol.u)
        end
        if save_end
            open(string("./",Extrastring,"Endstate_",sigmas[i],".serialjl"),"w") do io
                serialize(io,System.state)
            end
        end
    end
    writedlm(string("./",Extrastring,"Sigmaline.txt"),abs.(Orders))
    writedlm(string("./",Extrastring,"Phis.txt"),phis)
end
n = 1000#parse(Int64,ARGS[1])
startsigma = 1.0
if length(ARGS)>=2
    startsigma=parse(Float64,ARGS[2])
end
state = nothing
if length(ARGS)>=3
    open(ARGS[3]) do io
        global state = deserialize(io)
    end
end
Extrastring = "Up"
if length(ARGS)>=4
    Extrastring = ARGS[4]
end
Dothething(n,collect(startsigma:0.1:4.0),state=state,Extrastring=Extrastring)
