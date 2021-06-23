using DelimitedFiles, DifferentialEquations, Statistics, LinearAlgebra
#cd("/Users/simon3000/Documents/AG Schöll v2/simulation/FHN/MSF_FHN/testrun")

#network parameter
const sigma = parse(Float64,ARGS[1]) #MSF for adaptive diffusive networks
const H0=parse(Float64,ARGS[2]) #plasticity function
const DH0_e=parse(Float64,ARGS[3]) #plasticity function jacobian first entry

const epsi1 = 0.08 #LI07b
const epsi2 = 0.01 #coupling strength - membrane potential time separation
const iapp = 0 #LI07b
const a = 0.7 #LI07b
const b = 0.2 #LI07b
const B = 0 #LI07b
const om = 0 #LI07b
const v_shp = 0.05 #LI07b
const tau_syn = 1/1.2 #LI07b
const alpha0 = 2 #LI07b
const rho = 1 #adaption strength, MSF for adaptive diffusive networks
dynvar=4 #number of dynamical variables including coupling

println("values sigma=$(sigma), H0=$(H0), DH0=$(DH0_e)")

#plot parameter
xrange=4 #range x axis in px, factor 2!
yrange=4 #range y axis in px, factor 2!
xmid=1 #center x axis
ymid=0 #center y axis
res=100 #resolution, subdivision of integer steps

function save_log()
	writedlm("log_MSF_FHN_sigma$(sigma)_H0$(H0)_DH0$(DH0_e).csv",
	[
	dynvar "number of dynamical variables including coupling";
	H0 "H0 plasticity function";
	DH0_e "plasticity function first entry jacobian"
	xrange "range x axis in px";
	yrange "range y axis in px";
	xmid "center x axis";
	ymid "center y axis";
	res "resolution";
	sigma "sigma -> overall coupling strength";
	per "period of stable solution"
	], '\t')
end

function find_period(var1, var2, dt, transient) # finds average period in tunits, transient in timeunits
    phase = atan.(var2[Int(transient/dt):end], var1[Int(transient/dt):end]) # polar angle of complex number var1 + var2 * im
    period_length = Int64[] # vector for period lengths
    count = 0 # counts timesteps
    for i in eachindex(phase)[1:end-1]
        count += 1
        sign(phase[i]) != sign(phase[i+1]) && (push!(period_length, count); count = 0) # push counter if sign of phase(t) and phase(t+1) differ
    end
    deleteat!(period_length,1) # delete first incomplete period
	isodd(length(period_length)) && deleteat!(period_length,1) # delete another period if number of measured periods is odd (in case period halfs are not equal)
    avr_per = 2*sum(period_length)/length(period_length)
    return avr_per*dt
end

function FHN_model(du, u, t, p)
    #synchronized solution
    du[1] = 1/epsi1 * (u[1] - 1/3 * u[1]^3 - u[2] + iapp + rho * sigma * H0 * u[3] * u[1]) #membrane potential
    du[2] = u[1] + a - b * u[2] #recovery variable
    du[3] = alpha0 * (1 - u[3])/(epsi1 * (1+exp(-u[1]/v_shp))) - u[3]/tau_syn #synaptic variable

	#functions and jacobis
	Df11 = 1/epsi1*(1-u[1]^2)
    Df31 = alpha0*(1-u[3])*exp(-u[1]/v_shp)/(epsi1*v_shp*(1+exp(-u[1]/v_shp))^2)
    Df33 = -alpha0/(epsi1*(1+exp(-u[1]/v_shp)))-1/tau_syn
    Df = [Df11 -1/epsi1 0; 1 -b 0; Df31 0 Df33]
    DH0 = [DH0_e 0 0]
    Gs = [1/epsi1*u[1]*u[3], 0, 0]
    D1Gs = [1/epsi1*u[3] 0 0; 0 0 0; 0 0 0]
	D2Gs = [0 0 1/epsi1*u[1]; 0 0 0; 0 0 0]

    du[4:6] = Df*u[4:6]+sigma*rho*H0*(D1Gs*u[4:6]+D2Gs*(u[4:6]-mudr_r*u[4:6]+mudr_i*u[8:10]))-sigma*Gs*u[7] #real part zeta
    du[7] = -epsi2*((rho*DH0*(mudr_r*u[4:6]-mudr_i*u[8:10]))[1]+u[7]) #real part kappa

    du[8:10] = Df*u[8:10]+sigma*rho*H0*(D1Gs*u[8:10]+D2Gs*(u[8:10]-mudr_i*u[4:6]-mudr_r*u[8:10]))-sigma*Gs*u[11] #imaginary part zeta
    du[11] = -epsi2*((rho*DH0*(mudr_i*u[4:6]+mudr_r*u[8:10]))[1]+u[11]) #imaginary part kappa
end


mono_init = [if i+1 == j 1.0 else 0.0 end for i in 1:2*dynvar, j in 1:dynvar+1] #starting conditions, iterator i factor 2 bc of imag. parts of var. equation
floq_plane = Array{ComplexF64,2}(undef,xrange*res+1,yrange*res+1) #initialize floq matrix
sol_period = Array{Float64,1}[] #initialize array to save solutions

for j in 0:Int(xrange*res) #nested loop to construct mu, real
	global mudr_r=round(j/res-xrange/2+xmid,digits=6) #real part eigenvalue laplacian (Mu) Divided thru Rowsum
	for k in 0:Int(yrange*res) #imaginary
		global mudr_i=round(k/res-yrange/2+ymid,digits=6) #imaginary part eigenvalue laplacian (Mu) Divided thru Rowsum

		mono = Array{Complex, 1}[] #initialize monodromy matrix
		for i in 1:dynvar+1 #1st run find initial conditions on periodic orbit, 2nd-5th run computing monodromy matrix (4 dynamical variables)
			k==0 && j==0 || i==1 && continue #only calculate periodic orbit on very first run -> no changes for different mu

			# starting conditions
			if i==1
			    u0 = vcat([
			        0, # membrane potential sync
			        -0.5, # recovery variable sync
			        0.2], # synaptic variable
			        mono_init[:,i]) #initialize starting vector monodromy
			    tspan = (0.0, 10000.0)
			else u0 = vcat(start_sync, mono_init[:,i]); tspan=(0.0, per)
			end

			println("solving with: µ/r = $(complex(mudr_r,mudr_i)), var = $(i)")
			#println("u0=$(u0)")
			#println("tspan=$(tspan)")
			prob = ODEProblem(FHN_model, u0, tspan)
			#print("solving ODE (j=$(j), k=$(k), i=$(i))...")
			sol = solve(prob, Tsit5(), abtol=1e-6, reltol=1e-6, saveat=0.005)
			#println("length(sol) = $(size(sol))")
			if i == 1
				global start_sync = sol[1:3,end] # initial conditions on periodic orbit
				println("stable conditions initialized.")


				global per = find_period(sol[1,:],sol[2,:],0.005,9000)
				println("found periodic solution; T=$(per) t.units.")
				writedlm("sol_find_period_sigma$(sigma)_H0$(H0)_DH0$(DH0_e).csv", sol)
				writedlm("period_sigma$(sigma)_H0$(H0)_DH0$(DH0_e).csv", per)
			else
                push!(mono, complex.(sol[4:7,end],sol[8:11,end])) # save entry monodromy matrix
                push!(sol_period, vcat([j,k,i], sol[:,end])) # save solution after one period
			end
			#writedlm("sol_j$(j)_k$(k)_i$(i).csv", sol) #write in single file (sol_period) to save memory
		end

		floq = eigvals(hcat(mono...)) #floquet multiplicators
		#all(abs.(floq) .< 1) ? println("--STABLE--") : println("--UNSTABLE--")
		max_idx = findmax(abs.(floq))[2] # index of maximum absolute value
		floq_plane[k+1,j+1] = floq[max_idx] # save complex floquet multiplicator
	end
end
writedlm("floq_plane_sigma$(sigma)_H0$(H0)_DH0$(DH0_e).csv", floq_plane)
writedlm("sol_period_sigma$(sigma)_H0$(H0)_DH0$(DH0_e).csv", hcat(sol_period...))
save_log()
