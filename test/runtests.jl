using FirstOrderSolvers
using Base.Test

@testset "Tests:" begin
    @testset "testlbfgs.jl" begin
        include("testlbfgs.jl")
    end

    @testset "conjugateGradient.jl" begin
        include("conjugateGradient.jl")
    end

    @testset "HSDEAffine.jl" begin
        include("HSDEAffine.jl")
    end

    @testset "affinepluslinear.jl" begin
        include("affinepluslinear.jl")
    end

    @testset "testDRandGAPA.jl" begin
        include("testDRandGAPA.jl")
    end

    @testset "testsupermann.jl" begin
        include("testDRandGAPA.jl")
    end

    @testset "testPSD.jl" begin
        include("testPSD.jl")
    end

    @testset "testprint.jl" begin
        include("testprint.jl")
    end

    @testset "testconvex.jl" begin
        include("testconvex.jl")
    end
end
