using Test
include("KuramotoCustomIntegrator.jl")

@testset "Kuramoto-integration" begin
    @testset "Kuramoto-system" begin
        @test typeof(Kuramotosystem([1,2,3,4],1.,2.,3.,4.,5.,[1,2,3,4]))==Kuramotosystem
        @test length(ConstructStatefromPhis([1.,2.,3.,4.,5.,6.]))==(6+6*6)
        Testsystem = Kuramotosystem([1,2,3,4])
        @test Testsystem.state[1:4]==[1,2,3,4]
        @test typeof(Testsystem.state) == Array{Float64,1}
        @test Testsystem.α==0
        @test Testsystem.β==pi/2
        @test Testsystem.N==4
        @test Testsystem.σ==2.5
        @test Testsystem.ϵ==0.01
        @test Testsystem.omegas==zeros(Testsystem.N)
        Testsystem = Kuramotosystem(randn((50)),ϵ=1,β=pi)
        @test Testsystem.ϵ == 1
        @test -sin(Testsystem.state[4]-Testsystem.state[39]+Testsystem.β)==Testsystem.state[50*4+39]
    end
    @testset "Integrationfunctions" begin
        Testsystem = Kuramotosystem(randn(50))
        dstate = zeros(length(Testsystem.state))
        Derivativesphi(dstate,Testsystem)
        @test dstate[Testsystem.N+1:end] == zeros(Testsystem.N^2)
        dstate = zeros(length(Testsystem.state))
        Derivativeskappa(dstate,Testsystem)
        @test dstate[1:Testsystem.N] == zeros(Testsystem.N)
        Testsystem = Kuramotosystem(zeros(50))
        @test order(Testsystem)==1
        Testsystem = Kuramotosystem(ones(50),β=0,α=0)
        Derivativesphi(dstate,Testsystem)
        @test dstate[1:Testsystem.N] == zeros(Testsystem.N)
        Derivativeskappa(dstate,Testsystem)
        @test dstate[Testsystem.N+1:end] == zeros(Testsystem.N^2)
        Testsystem = Kuramotosystem(LinRange(0,2pi,51)[1:end-1])
        @test order(Testsystem)+1≈1
        TestPhis = zeros((50,100))
        @test order(TestPhis)==ones(100)
        @test order(TestPhis[1:25,74:90])==ones(17)
        @test typeof(CustomRK4(Testsystem,(0.,1.)))==NamedTuple{(:t,:u),Tuple{Array{Float64,1},Array{Float64,2}}}
        @test_throws ArgumentError CustomRK4(Testsystem,(0.,10.),save_indices="gfa")
        @test_throws BoundsError CustomRK4(Testsystem,(0.,10.),save_indices=length(Testsystem.state)+1)
        @test_throws MethodError CustomRK4(Testsystem,(0.,10.),save_times=0.33)
        @test_throws ErrorException CustomRK4(Testsystem,(0.,10.),save_times=[0.33])
        test_stuff = CustomRK4(Testsystem,(1.,11.),save_times = [8], save_indices = 2)
        @test test_stuff.t == [1,8]
        @test size(test_stuff.u)==(1,2)
    end
end
