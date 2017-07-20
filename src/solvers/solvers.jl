export GAP, GAPA, GAPP,  DR, AP, Dykstra, FISTA

include("gap.jl")
include("gapa.jl")
include("gapproj.jl")
include("dykstra.jl")
include("fista.jl")

#Some shortnames for common algorithms
DR(α=0.5; kwargs...) = GAP(α, 2.0, 2.0; kwargs...)
AP(α=1; kwargs...)   = GAP(α, 1.0, 1.0; kwargs...)
