using Plots
using DelimitedFiles
using Random
include("KuramotoCustomIntegrator.jl")


function SpecialforthisspecialtestCustomRK4(system, tspan; dt = 0.1, save_indices=1:system.N, save_times=collect(tspan[1]:dt:tspan[2]))
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
    save_times_iterator = 1
    if tspan[1]∈save_times
        save_times_iterator +=1
    end
    if save_times_iterator==1
        saved_stuff = (t = [tspan[1], save_times...],
        u = zeros(length(system.state[save_indices]), length([tspan[1], save_times...])),)
    elseif save_times_iterator==2
        saved_stuff = (t = save_times,
        u = zeros(length(system.state[save_indices]), length(save_times)),)
    end
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
        MeanKappas = zeros(2)
        for i in 1:2
            for k in (collect(1:319),collect(320:400))[i]
                for l in (collect(1:319),collect(320:400))[i]
                    MeanKappas[i]+=system.state[system.N*k+l]
                end
            end
            MeanKappas[i]/=length((collect(1:319),collect(320:400))[i])*length((collect(1:319),collect(320:400))[i])
        end
        for i in 1:2
            for k in (collect(1:319),collect(320:400))[i]
                for l in (collect(1:319),collect(320:400))[i]
                    system.state[system.N*k+l]=MeanKappas[i]
                end
            end
        end
        if integrated_times[i]==save_times[save_times_iterator]
            saved_stuff.u[:, save_times_iterator] .= system.state[save_indices]
            if save_times_iterator < length(save_times)
                save_times_iterator+=1
            end
        end
        if rem(100 * (i-1), length(tspan[1]:dt:tspan[2])-1) == 0
            println(100 * (i-1) / (length(tspan[1]:dt:tspan[2])-1))
        end
    end
    return saved_stuff
end


Oszis = 400
cutoff = 0.8
s = (3.2)
t = (0.,3000.)

Random.seed!(314159265)
omegas = sort!(2*rand(Oszis)-1. .*ones(Oszis))
Phis = zeros(Oszis)
for i in round(Int, cutoff*length(Phis), RoundDown):length(Phis)
    Phis[i]+=0.3*pi
end
System = Kuramotosystem(Phis,omegas=omegas,α=0.15*pi, β=-0.53*pi, σ = s)
for i in round(Int, cutoff*length(Phis), RoundDown):length(Phis)
    for j in 1:round(Int, cutoff*length(Phis),RoundDown)-1
        System.state[Oszis+Oszis*(i-1)+j]=0
        System.state[Oszis+Oszis*(j-1)+i]=0
    end
end
System.state[Oszis+1:end].+= 0.01*rand.(Oszis*Oszis)
full_sol = CustomRK4(System, t, save_indices=1:Oszis+Oszis^2, dt = 0.1, save_times=collect(2990:0.5:3000))
#open(string("./Julia/FullSims/",s,"avgkappa.txt"),"w")  do io
#    writedlm(io,full_sol)
#end
