using LightGraphs, Plots, ColorSchemes
using DifferentialEquations
using Plots
using DelimitedFiles, FFTW, Random
using Serialization
using Statistics
pyplot()

function CreateOrders(pathtofiles;n = collect(1:50),sigmas = collect(0.0:0.01:10.0))
    Orders = zeros((length(n),length(sigmas)))
    for i in 1:length(n)
        for j in 1:length(sigmas)
            synchmatrix = readdlm(string(pathtofiles,sigmas[j],"_",n[i],"_ordermatrix.csv"))
            Orders[i,j] = mean(synchmatrix)
        end
    end
    open(string(pathtofiles,"Orders.txt"),"w") do io
        writedlm(io,Orders)
    end
end

function CreateCliques(pathtofiles;n = collect(1:50),sigmas = collect(0.0:0.01:10.0))
    Cliques = Array{Array{Any}}(undef,length(n),length(sigmas))
    for i in 1:length(n)
        for j in 1:length(sigmas)
            Orders = readdlm(string(pathtofiles,sigmas[j],"_",n[i],"_ordermatrix.csv"))
            g = SimpleGraph(length(n))
            for k in 1:length(n)
                for l in 1:length(n)
                    if Orders[k,l]==1
                        add_edge!(g, k, l)
                    end
                end
            end
            Cliques[i,j]=maximal_cliques(g)
        end
    end
    serialize(string(pathtofiles, "CliquesSerialized.serialjl"),Cliques)
    open(string(pathtofiles,"Cliques.txt"),"w") do io
        writedlm(io,Cliques)
    end
end
function CreatePic(pathtofiles,pathtopicture;n = collect(1:50),sigmas = collect(0.0:0.01:10.0),picxticks = (0:100:1000,sigmas[1:100:1001]),pictitle = nothing)
    pyplot()
    Cliques = deserialize(string(pathtofiles, "CliquesSerialized.serialjl"))
    Orders = readdlm(string(pathtofiles, "Orders.txt"))
    heatmap(Orders,ylabel="Oscillators in Cluster 1 at start", xlabel="Sigma", xticks=picxticks,title=pictitle)
    Expected = zeros(size(Cliques))
    for i in 1:length(n)
        for j in 1:length(sigmas)
            if n[i]<50
                Expected[i,j] = length(Cliques[i,j])==2 && (length(Cliques[i,j][1])==n[i] || length(Cliques[i,j][2])==n[i])
            end
            if n[i]==50
                Expected[i,j] = length(Cliques[i,j])==1 && length(Cliques[i,j][1])==n[i]
            end
        end
    end
    contour!(Expected, legend = false, color="blue", linewidth=0.1)
    savefig(pathtopicture)
end
function ClusterScatterPlot(Phis,n::Integer)
    scatter(cos.(Phis[1:n]),sin.(Phis[1:n]),c=:yellow,markerstrokewidth=0)
    scatter!(cos.(Phis[n+1:end]),sin.(Phis[n+1:end]),c=:blue,markerstrokewidth=0)
    plot!(xlims = [-1,1],ylims=[-1,1],legend=false)
end
function CreateMovieOfOscis(Phis,n,pathtomovie;fps=100)
    Size = size(Phis)
    anim = Animation()
    for i in 1:Size[2]
        ClusterScatterPlot(Phis[:,i],n)
        frame(anim)
    end
    gif(anim,pathtomovie,fps=fps)
end
function CreateFFTDiagram(Orders;error=true)
    if error
        F = abs.(fft(Orders[2:end-1].-mean(Orders[2:end-1]))) |> fftshift
    else
        F = abs.(fft(Orders.-mean(Orders))) |> fftshift
    end
    freq = fftfreq(length(F),10) |> fftshift
    plot(freq,F,xlim=[-1,1]),F,freq
end
function OrderFFT(runname,n,sigma;error=true,prefix="Up")
    Orders = readdlm(string("./Data/",runname,"/Data/",n,"/",prefix,"Timeline_",sigma,".txt"),'\t',Complex{Float64})
    CreateFFTDiagram(Orders,error=error)
end
function RemoveJumps(x,y,cutoff)
    x1 = copy(x)
    y1 = copy(y)
    Indices = FindJumps(y1,cutoff)
    for i in Indices
        insert!(x1,i,NaN)
        insert!(y1,i,NaN)
    end
    x1,y1
end
function FindJumps(y,cutoff;multiplier = 1.5,upperlimit=0.15)
    Indices = Array{Int64}(undef,0)
    for i in 3:length(y)
        if y[i]>y[i-1] && (abs(y[i]-y[i-1])>upperlimit || (abs(y[i]-y[i-1])>cutoff && abs(y[i]-y[i-1])>abs(y[i-1]-y[i-2])*multiplier))
            push!(Indices,i)
        end
    end
    Indices
end
function PlotLineWithJumps(x,y,Indices,c;label="",alpha=0.5,kwargs...)
    if Indices[1] != 1
        insert!(Indices,1,1)
    end
    for index in 1:length(Indices)-1
        plot!(x[Indices[index]:Indices[index+1]-1],y[Indices[index]:Indices[index+1]-1],c=c,label="";kwargs...)
        plot!(x[Indices[index+1]-1:Indices[index+1]],y[Indices[index+1]-1:Indices[index+1]],alpha=alpha,c=c,label="",style = :dot;kwargs...)
    end
    plot!(x[Indices[end]:end],y[Indices[end]:end],c=c,label=label;kwargs...)
end
function PlotBifDiagram2(sigmas,Orders;alpha = 0.5,cutoff=0.01,multiplier=1.5,kwargs...)
    shape = size(Orders)
    for i in 1:shape[1]
        Indices = FindJumps(Orders[i,1:end],cutoff,multiplier=1.5)
        PlotLineWithJumps(sigmas,Orders[i,1:end],Indices,get(ColorSchemes.bamako,i,(1,shape[1]));alpha = alpha,kwargs...)
    end
    plot!()
end
function PlotBifDiagram(sigmas,Orders;skips=1,downsweep=true,cutoff=0.01,kwargs...)
    n=size(Orders)[1]
    if downsweep
        n = Int(size(Orders)[1]/2)
    end
    for i in 1:n
        previndex = 1
        for j in (skips+1):length(sigmas)
            if abs(Orders[i,j+skips]-Orders[i,j+skips-1])>cutoff
                plot!(sigmas[previndex:j-1],Orders[i,previndex+skips:j+skips-1],c=get(ColorSchemes.bamako,Orders[i,1],(0.5,1.0)),label="";kwargs...)
                previndex = j
            end
        end
        plot!(sigmas[previndex:end],Orders[i,previndex+skips:end],c=get(ColorSchemes.bamako,Orders[i,1],(0.5,1.0)),label=string("n_1 = ",Orders[i,1]);kwargs...)
        previndex = 1
        if downsweep
            plot!(sigmas[1:end],Orders[i+n,skips+1:end],c=get(ColorSchemes.bamako,Orders[i+n,1],(0.5,1.0)),label="",style=:dash;kwargs...)
        end
    end
    plot!(xlabel = "Sigma", ylabel="2nd Orderparameter")
end
function PhisvsPhis(OszisPhis,CollPhis,n)
    scatter(OszisPhis[1:n],CollPhis[1:n],markerstrokewidth=0,c=:blue,msize=2)
    scatter!(OszisPhis[n+1:end],CollPhis[n+1:end],markerstrokewidth=0,c=:green,msize=2.0)
    Linecoords = [round(OszisPhis[1]*10-1)/10,round(OszisPhis[end]*10)/10]
    plot!(Linecoords,Linecoords,c=:black,linewidth=0.2)
    plot!(legend=false,xlabel="\$\\phi_i\$",ylabel="\$\\^\\phi_i\$")
end
function Kappadiffs(OszisKappas,CollKappas,n;N=1000)
    OszisKappas=reshape(OszisKappas,(N,N))
    CollKappas=reshape(CollKappas,(2,2))
    Kappadiffs=zeros((N,N))
    for i in 1:N
        for j in 1:N
            Kappadiffs[i,j]=OszisKappas[i,j]-CollKappas[i<=n ? 1 : 2,j<=n ? 1 : 2]
        end
    end
    p=heatmap(Kappadiffs)
    return p,Kappadiffs
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
function SortState(state,params,perm1)
    k = reshape(state[params.N+1:end],(params.N,params.N))[perm1,perm1]
    [state[1:params.N][perm1]; [k...]]
end
function CalcPerm(state,params)
    perm1 = sortperm(mod2pi.(state[1:params.N]))
    state = SortState(state,params,perm1)
    dstate = zeros(params.N)
    params.omegas[1:end].=params.omegas[p]
    Derivativesphi2(dstate,state,params,0)
    perm1 = sortperm(dstate)
    SortState(state,params,perm1)
end
function CreateParams2(omegas, indices; α=0*pi, β=-0.53*pi, σ=3.5, ϵ=0.01)
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
function PlotHeatmap(state,seed,sigma)
    heatmap(reshape(state[51:end],(50,50)),xlabel="\$i\$",ylabel="\$j\$",color=:bwr,cbtitle="\$\\kappa_{ij}\$",clims=(-1,1))
    savefig(string("./Pics/FirstPicture/Coupling_",seed,"_",sigma,".png"))
    savefig(string("./Pics/FirstPicture/Coupling_",seed,"_",sigma,".svg"))
    heatmap(reshape(abs.(state[51:end]),(50,50)),xlabel="\$i\$",ylabel="\$j\$",color=:bwr,cbtitle="\$\\left|\\kappa_{ij}\\right|\$",clims=(0,1))
    savefig(string("./Pics/FirstPicture/Coupling_",seed,"_",sigma,"_abs.png"))
    savefig(string("./Pics/FirstPicture/Coupling_",seed,"_",sigma,"_abs.svg"))
end
function CreateHeatmap(seed,sigma)
    Random.seed!(seed)
    omegas = sort(0.5 .* rand(50) .- 0.25)
    params = CreateParams2(omegas,[collect(1:50),],σ=sigma)
    state = deserialize("./Data/50OszisSweep/Data/-53_$(seed)_$(sigma)_state.serialjl")
    state = CalcPerm(state,params)
    PlotHeatmap(state,seed,sigma)
end
function StabilitySweep(stability)
    x = collect(0:0.01:5)
    y = collect(1:50)
    heatmap(x,y,stability,ylabel="Oscillators in Cluster 1",xlabel="Sigma",cbar_title="LLE")
end
function Frequencyhistogram(freqs;nbins=100, left_bin_edges=range(-0.25,0.25,length=nbins+1)[1:end-1])
    normalization,length = size(freqs)
    hist = zeros(nbins,length)
    for i in 1:length
        for value in freqs[1:end,i]
            bin = findlast(x->value>x,left_bin_edges)
            hist[bin,i]+=1
        end
    end
    hist
end
function Frequencies(freqs;nbins=100, left_bin_edges=range(-0.25,0.25,length=nbins+1)[1:end-1], transformation = identity)
    normalization,length = size(freqs)
    hist = Frequencyhistogram(freqs,nbins=nbins,left_bin_edges=left_bin_edges)
    heatmap(collect(0:0.01:10),left_bin_edges,transformation(hist ./normalization),cbar_title="Probability",xlabel="sigma",ylabel="frequency")
end
function Vockpicture(freqs,n)
    for i in 1:n
        freqs[50*(i-1)+1:50*(i),1:end] .-= mean(freqs[50*(i-1)+1:50*i,end])
    end
end
function Orderparameter(frequencies)
    Matrix = zeros(50,50)
    for i in 1:50
        for j in 1:50
            Matrix[i,j]=abs(frequencies[i]-frequencies[j])<0.001
        end
    end
    Matrix
end
function CountFrequencies(freqs)
    count = [freqs[1]]
    for f in freqs
        bla = true
        for c in count
            if abs(c-f)<0.001
                bla = false
            end
        end
        if bla
            count = [count f]
        end
    end
    count
end
function GetOrders(filepath)
    dir = readdir(filepath,join=true)
    counter = 0
    for file in dir
        if file[end-5:end]=="rs.txt"
            counter +=1
        end
    end
    Orders = zeros(counter,1001)
    counter = 1
    for file in dir
        if file[end-5:end]=="rs.txt"
            bla = readdlm(file)
            Orders[counter,1:end] = bla
            counter +=1
        end
    end
    Orders
end
function plotfigure1(order1,order2,order3)
    freqs = readdlm("./Data/Eckdaten/Frequencymatrix.txt")
    Pointer = readdlm("./Data/Eckdaten/2C_1C_whatever_Pointer.txt")
    keys = sortperm(sum(order1[1:end,1:600],dims=2)[1:end])
    Random.seed!(0)
    omegas1 = sort(0.5*rand(50) .- 0.25)
    Random.seed!(18)
    omegas2 = sort(0.5*rand(50) .- 0.25)
    plot()
    p1 = PlotBifDiagram2(collect(0:0.01:10),order1[keys,1:end];xlabel="\$\\sigma\$",xticks=[0,6],yticks=[0,1],ylabel = "\$\\overline{R}_2\$",legend=false,xlim=(0,6))
    annotate!(0.1,1,Plots.text("a",:top,:left,8))
    annotate!(3.5,minimum(order1[findall(isequal(1.0),Pointer[1:30]),350:end])-0.01,Plots.text("gradual",:top,:left,8))
    annotate!(3,minimum(order1[1:end,300:end])-0.01,Plots.text("explosive",:top,:left,8))
    plot()
    p2 = PlotBifDiagram2(0:0.01:10,order2;xlabel="\$\\sigma\$",xticks=[0,6],yticks=[0,1],ylabel = "\$\\overline{R}_2\$",legend=false,xlim=(0,6))
    annotate!(0.1,1,Plots.text("b",:top,:left,8))
    plot()
    p3 = PlotBifDiagram2(0:0.01:10,order3;xlabel="\$\\sigma\$",xticks=[0,6],yticks=([0,1],["",""]),legend=false,xlim=(0,6))
    annotate!(0.1,1,Plots.text("c",:top,:left,8))
    scatter!(omegas1,c=get(ColorSchemes.bamako,0),xlabel="i",ylabel="\$\\omega_i\$",guidefontsize=1,annotations=(1,0.2,Plots.text("d",:top,:left,8)),ylims=(-0.25,0.25),yticks=([-0.25,0.25],["\$-\\^{\\omega}\$","\$\\^{\\omega}\$"]),xticks=[1,50],markerstrokewidth=0,markersize=1.5,label="",inset=(1,bbox(0.4,0.35,0.5,0.5)),subplot=2)
    scatter!(omegas2,c=get(ColorSchemes.bamako,0.8),ylims=(-0.25,0.25),label="",subplot=2,markerstrokewidth=0,markersize=1.5)
    plot()
    #p4 = scatter(omegas1,c=get(ColorSchemes.bamako,0),xlabel="Index i",ylabel="\$\\omega_i\$",annotations=(1,0.2,Plots.text("d",:top,:left,8)),ylims=(-0.25,0.25),yticks=[-0.25,0.25],xticks=[1,50],label="")
    #scatter!(omegas2,c=get(ColorSchemes.bamako,1),label="c)",legend=:bottomright)
    l = @layout [a{0.5h};b c;d e]#;_ d{0.5w} _]
    l = @layout [a{0.6h};b c]
    #p_all = plot(p1,p2,p3,p4,plot([0],[0]),plot([0],[0]),layout=l,size=(600,800))
    p_all = plot(p1,p2,p3,layout=l,size=(350,350),link=:none,guidefontsize=10,grid=false,tickfontsize=8)
    plot!(p3,subplot=2,guidefontsize=8,tickfontsize=7)
    p_all
end
function Figure1Order()
    order1 = GetOrders("./Data/50OszisKuramotos/Data/")
    order2 = GetOrders("./Data/50OszisKuramotos/Data/0/")
    order3 = GetOrders("./Data/50OszisKuramotos/Data/18/")
    plotfigure1(order1,order2,order3)
end
function Figure1Synch()
    order1 = transpose(readdlm("./Data/50OszisSweep/Orders.txt"))
    order2 = GetOrders("./Data/50OszisSweepFixedSeed/Data/0/")
    order3 = GetOrders("./Data/50OszisSweepFixedSeed/Data/18/")
    plotfigure1(order1,order2,order3)
end
function Figure3()
    sigmas = collect(0:0.01:5)
    sigmas2 = collect(0:0.1:4)
    ns = collect(500:50:1000)
    Orders = readdlm("./Data/Singleruns/Orders.txt")
    Orders2 = readdlm("./Data/FullSystemOrders/Orders.txt")
    Orders[1:10,2:52]=Orders[1:10,52:-1:2]
    plot()
    plot!([100/32/sin(0.53*pi),100/32/sin(0.53*pi)],[0,1],c=:black,alpha=0.5, label="")
    plot!([100/32/sin(0.53*pi),5],[1,1],legend=false,alpha=0,fill=0,fillalpha=0.3,fillcolor=:black,fillstyle=:/)
    PlotBifDiagram(sigmas,Orders,cutoff=1;style=:dash)
    PlotBifDiagram(sigmas2,Orders2;downsweep=false, cutoff = 1, label="")
    plot!([100/32/sin(0.53*pi),100/32/sin(0.53*pi)],[0,1],c=:black,alpha=0.5, label="")
    plot!(xlabel="\$\\sigma\$",ylabel="\$\\overline{R}_2\$")
    CollPhis = deserialize("./Data/Comparisons/Collective_Phis.serialjl")
    OszisPhis = mod2pi.(deserialize("./Data/Comparisons/UpEndstate_2.0.serialjl")[1:1000])
    CollPhis = CollPhis .-CollPhis[1].+OszisPhis[1]
    n=500

    scatter!(OszisPhis[1:n],CollPhis[1:n],markerstrokewidth=0,c=:blue,msize=2.0,inset=(1,bbox(0.2,0.5,0.4,0.4)),subplot=2,xticks=:none,yticks=:none,legend=false)
    scatter!(OszisPhis[n+1:1000],CollPhis[n+1:end],markerstrokewidth=0,c=:green,msize=2.0,subplot=2,xticks=:none,yticks=:none,legend=false)
    Linecoords = [round(OszisPhis[1]*100-1)/100,round(OszisPhis[end]*100+1)/100]
    plot!(Linecoords,Linecoords,c=:black,linewidth=0.2,subplot=2)
    plot!(legend=false,xlabel="\$\\phi_i\$",ylabel="\$\\^\\phi_i\$",subplot=2)
    #plot(p2,inset=(p1,bbox(0.2,0.5,0.45,0.45)))
    #p1
end


function PlotHeatmapFromFile(filestring;kwargs...)
    state = deserialize(filestring)
    heatmap(reshape(state[51:end],(50,50)),c=:bwr,clims=(-1,1);kwargs...)
end
Pointer = readdlm("./Data/Eckdaten/2C_1C_whatever_Pointer.txt")
Synchs = readdlm("./Data/Eckdaten/Synchronizationparameter_101runs_50Oszis.txt")
blank1 = quiver([0.06],[0],quiver=([0.9],[0]),c=:black,framestyle=:none,xlim=(0,1),ylim=(-0.5,0.5),annotations=(0.05,0,Plots.text("\$\\sigma\$",:center,:right)))
blank2 = plot(foreground_color_subplot=:white)
l = @layout [a{0.08h};Plots.grid(2,4) b{0.075w}]
p_all = scatter!([0,0],[0,0],zcolor=[clims[1],clims[2]],label="",markerstrokecolor=:transparent,clims=clims,background_color_subplot=:transparent,framestyle=:none,inset=bbox(0.1,0,0.3,0.9,:center,:right),subplot=10)

#=
Pointer = zeros(101)
for i in 0:100
    frequencies = readdlm(string(path,"-53_$(i)_Frequencies.txt"))
    j=findlast(x->x<1,Synchs[i+1,1:end])
    Clusterf = CountFrequencies(frequencies[1:end,j])
    if Synchs[i+1,j]>0.75
        Pointer[i+1]=1.
    elseif length(Clusterf)==2
        Pointer[i+1]=0.
    else
        Pointer[i+1]=2
    end
end
Firstindex = collect(0:100)[Pointer .==2][1]
freqlow = readdlm("./Data/50OszisFrequencies_randomseed/Data/-53_$(Firstindex)_Frequencies.txt")
for i in collect(0:100)[Pointer .== 2][2:end]
    global freqlow = [freqlow; readdlm("./Data/50OszisFrequencies_randomseed/Data/-53_$(i)_Frequencies.txt")]
end
=#
#Folders = ["000+053","-030-053","000-053","030+053","-030+053","030-000","-030-000"]
#for loc in Folders
#    CreatePic(string("./Data/",loc,"/Data/"),string("./Pics/Phillipescans/",loc),pictitle=loc)
#end
#==
CreateCliques("./Data/000-053_0/Data/")
CreateOrders("./Data/000-053_0/Data/")
CreatePic("./Data/000-053_0/Data/","./Pics/Phillipescans/Omega0")
==#
