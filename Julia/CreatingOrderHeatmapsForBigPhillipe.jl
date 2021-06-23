using LightGraphs
using DifferentialEquations
using Plots
using DelimitedFiles
using Serialization
using Statistics

pyplot()
n = collect(1:50)
sigma = collect(0.:0.01:10.)
#
# Heatmaps with the cliques
Cliques = deserialize("./Data/BigPhillipeScan/CliquesSerialized.serialjl")
Orders = readdlm("./Data/BigPhillipeScan/Orders.txt")
heatmap(transpose(Orders),ylabel="Oscillators in Cluster 1 at start", xlabel="Sigma", xticks=(0:100:1000,sigma[1:100:1001]))
#heatmap(map(length,Cliques),title="Number of synchronized clusters",ylabel="Oscillators in Cluster 1 at start", xlabel="Sigma", xticks=(0:100:1000,sigma[1:100:1001]))
#plot!([0,1000],[25,25],color="red",legend=false)
#savefig("RandomW_NumberofClusters")
#heatmap(map(x->x>5 ? 0 : x,map(length,Cliques)),title="Number of synchronized clusters, >5->0",ylabel="Oscillators in Cluster 1 at start", xlabel="Sigma", xticks=(0:100:1000,sigma[1:100:1001]))
#plot!([0,1000],[25,25],color="red",legend=false)
#savefig("RandomW_NumberofClustersTruncated")
Expected = zeros(size(Cliques))
for i in 1:length(n)
    for j in 1:length(sigma)
        if n[i]<50
            Expected[i,j] = length(Cliques[i,j])==2 && (length(Cliques[i,j][1])==n[i] || length(Cliques[i,j][2])==n[i])
        end
        if n[i]==50
            Expected[i,j] = length(Cliques[i,j])==1 && length(Cliques[i,j][1])==n[i]
        end
    end
end
#heatmap(Expected,title="Same end as beginning?",ylabel="Oscillators in Cluster 1 at start", xlabel="Sigma", xticks=(0:100:1000,sigma[1:100:1001]))
#plot!([0,1000],[25,25],color="red",legend=false)
contour!(Expected, legend = false, color="blue", linewidth=0.1)
savefig("RandomW_OverlayedOrders")
#=
orders = readdlm("./Data/BigPhillipeScan/Orders.txt")
scatter(sigma[101:201],orders[126:201,1:end],legend=false,xlabel="Sigma",ylabel="Orderparameter")
savefig("RandomW_ScatterZoom")
plot(sigma[101:201],orders[126:201,1:end],legend=false,xlabel="Sigma",ylabel="Orderparameter")
savefig("RandomW_LineZoom")
=#
