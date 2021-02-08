using DelimitedFiles
using Statistics
include("KuramotoCustomIntegrator.jl")

Sigmas = collect(2.5:0.1:4.0)
indices = (collect(1:319),collect(320:400))
for s in Sigmas
    Data = readdlm(string("./Julia/FullSims/",s,".txt"),Float64,skipstart=1)
    Data = reshape(Data,(160400,402))
    Omegas = zeros(400)
    for i in 1:length(Omegas)
        Omegas[i]=(Data[i,end-1]-Data[i,1])/100
    end
    open(string("./Data/Omegas_",s,".txt"),"w")  do io
        writedlm(io,Omegas)
    end
    Orders = [zeros(401),zeros(401)]
    Orders[1] = order(Data[1:319,1:end-1])
    Orders[2] = order(Data[320:end,1:end-1])
    open(string("./Data/Orders_",s,".txt"),"w")  do io
        writedlm(io,Orders)
    end
    Kappas = [zeros(401) for i in 1:4]
    for t in 1:401
        for i in 1:length(indices)
            for j in 1:length(indices)
                for k in indices[i]
                    for l in indices[j]
                        Kappas[2*i-2+j][t]+=Data[k+400*l,t]
                    end
                end
                Kappas[2*i-2+j][t]/=length(indices[i])*length(indices[j])
            end
        end
    end
    open(string("./Data/Kappas_",s,".txt"),"w") do io
        writedlm(io,Kappas)
    end
end
