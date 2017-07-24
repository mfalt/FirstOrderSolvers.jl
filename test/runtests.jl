using FirstOrderSolvers
using Base.Test

info("Test: conjugateGradient.jl")
include("conjugateGradient.jl")

info("Test: HSDEAffine.jl")
include("HSDEAffine.jl")

info("Test: affinepluslinear.jl")
include("affinepluslinear.jl")

info("Test: testDRandGAPA.jl")
include("testDRandGAPA.jl")

info("Test: testPSD.jl")
include("testPSD.jl")

info("Test: testprint.jl")
include("testprint.jl")

info("Test: testconvex.jl")
include("testconvex.jl")
