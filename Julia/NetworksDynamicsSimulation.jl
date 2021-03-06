using LightGraphs
using DifferentialEquations
using NetworkDynamics
using Plots#, LaTeXStrings

N = 100 # number of nodes
g = complete_graph(N)

#=
  Berner, Rico, Eckehard Schöll, and Serhiy Yanchuk.
  "Multiclusters in Networks of Adaptively Coupled Phase Oscillators."
  SIAM Journal on Applied Dynamical Systems 18.4 (2019): 2227-2266.
=#


@inline function kuramoto_plastic_edge!(de, e, v_s, v_d, p, t)
    # Source to Destination coupling

    # The coupling function is modeled by a differential algebraic equation with mass matrix 0
    # 0 * de[1] = e[2] * sin(v_s[1] - v_d[1] + α) / N - e[1] is equivalent to e[1] = e[2] * sin(v_s[1] - v_d[1] + α) / N

    de[1] =  e[2] * sin(v_s[1] - v_d[1] + α) / N - e[1]
    de[2] = - ϵ * (sin(v_s[1] - v_d[1] + β) + e[2])
    # Destination to source coupling
    # since the coupling function is not symmetric we have to compute the other direction as well

    de[3] =  e[4] * sin(v_d[1] - v_s[1] + α) / N - e[3]
    de[4] = - ϵ * (sin(v_d[1] - v_s[1] + β) + e[4])
    nothing
end

@inline function kuramoto_plastic_vertex!(dv, v, e_s, e_d, p, t)
    dv[1] = 0
    for e in e_s
        dv[1] -= e[1]
    end
    for e in e_d # other direction is stored at other index
        dv[1] -= e[3]
    end
end

saved_values      = SavedValues(Float64, NTuple{N,Float64})
cb = SavingCallback((u,t,integrator)->(Tuple(u[1:N])), saved_values, saveat=180:200)
    # Parameter definitions
ϵ = 0.1
α = .2π
β = -.95π

# NetworkDynamics Setup
plasticvertex = ODEVertex(f! = kuramoto_plastic_vertex!, dim =1)
mass_matrix_plasticedge = zeros(4,4)
mass_matrix_plasticedge[2,2] = 1. # 1st and 3rd internal varibale are set to 0
mass_matrix_plasticedge[4,4] = 1.

plasticedge = ODEEdge(f! = kuramoto_plastic_edge!, dim=4, sym=[:es, :ks,:ed,:kd], mass_matrix = mass_matrix_plasticedge);
kuramoto_plastic! = network_dynamics(plasticvertex, plasticedge, g)

# ODE Setup & Solution
x0_plastic        = vcat(randn(N), ones(4ne(g)))
tspan_plastic     = (0., 200.)
params_plastic    = (nothing, nothing)
prob_plastic      = ODEProblem(kuramoto_plastic!, x0_plastic, tspan_plastic, params_plastic)
#=sol_plastic       = =#solve(prob_plastic, Rosenbrock23(), abstol = 1e-3, reltol = 1e-3, callback=cb, save_everystep=false, save_start=false, save_end = false)
println("I am done")
# Plotting
#v_idx = idx_containing(kuramoto_plastic!, :v)
plot_values = [[saved_values.saveval[j][i] for j in 1:length(saved_values.saveval)] for i in 1:length(saved_values.saveval[1])]
plot(saved_values.t,plot_values)

# Shows Coupling terms
#e_idx = idx_containing(kuramoto_plastic!, :k)
# plot!(sol_plastic, vars=e_idx, legend=false, color=:black, linewidth=0.001)
