using FirstOrderSolvers
using Base.Test

@testset "Tests:" begin
    info("Test: testlbfgs.jl")
    @testset "testlbfgs.jl" begin
        include("testlbfgs.jl")
    end
    info("Test: conjugateGradient.jl")
    @testset "conjugateGradient.jl" begin
        include("conjugateGradient.jl")
    end
    info("Test: HSDEAffine.jl")
    @testset "HSDEAffine.jl" begin
        include("HSDEAffine.jl")
    end
    info("Test: affinepluslinear.jl")
    @testset "affinepluslinear.jl" begin
        include("affinepluslinear.jl")
    end
    info("Test: testDRandGAPA.jl")
    @testset "testDRandGAPA.jl" begin
        include("testDRandGAPA.jl")
    end
    info("Test: testsupermann.jl")
    @testset "testsupermann.jl" begin
        include("testsupermann.jl")
    end
    info("Test: testPSD.jl")
    @testset "testPSD.jl" begin
        include("testPSD.jl")
    end
    info("Test: testprint.jl")
    @testset "testprint.jl" begin
        include("testprint.jl")
    end
    info("Test: testconvex.jl")
    @testset "testconvex.jl" begin
        include("testconvex.jl")
    end
end
