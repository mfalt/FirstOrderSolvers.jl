using FirstOrderSolvers
using LinearAlgebra
using Test
using Random

println("Test: conjugateGradient.jl")
include("conjugateGradient.jl")

println("Test: HSDEAffine.jl")
include("HSDEAffine.jl")

println("Test: affinepluslinear.jl")
include("affinepluslinear.jl")

println("Test: testDRandGAPA.jl")
include("testDRandGAPA.jl")

println("Test: testPSD.jl")
include("testPSD.jl")

println("Test: testprint.jl")
include("testprint.jl")

#println("Test: testspecific.jl")
#include("testspecific.jl")

# println("Test: testconvex.jl")
# include("testconvex.jl")
