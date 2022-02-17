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
