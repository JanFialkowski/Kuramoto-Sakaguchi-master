using LightGraphs
using DifferentialEquations
using Plots
using DelimitedFiles
using Serialization
using Statistics

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
#Folders = ["000+053","-030-053","000-053","030+053","-030+053","030-000","-030-000"]
#for loc in Folders
#    CreatePic(string("./Data/",loc,"/Data/"),string("./Pics/Phillipescans/",loc),pictitle=loc)
#end
CreateCliques("./Data/000-053_0/Data/")
CreateOrders("./Data/000-053_0/Data/")
CreatePic("./Data/000-053_0/Data/","./Pics/Phillipescans/Omega0")
