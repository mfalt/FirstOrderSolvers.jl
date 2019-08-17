"""
`GAPA(α=1.0, β=0.0; kwargs...)``

GAP Adaptive
Same as GAP but with adaptive `α1,α2` set to optimal `αopt=2/(1+sinθ)`
according to the estimate of the Friedrischs angle `θ`.
`β` is averaging factor between `αopt` and `2`: `α1=α2= (1-β)*αopt + β*2.0`.
"""
mutable struct GAPA  <: FOSAlgorithm
    α::Float64
    β::Float64
    direct::Bool
    options
end
GAPA(α=1.0, β=0.0; direct=false, kwargs...) = GAPA(α, β, direct, kwargs)

mutable struct GAPAData{T1,T2} <: FOSSolverData
    α12::Float64
    tmp1::Array{Float64,1}
    tmp2::Array{Float64,1}
    constinit::Base.RefValue{Bool} # If constant term has been calculated (linesearch feature)
    S1::T1
    S2::T2
end
#GAPAData{T1,T2}(α12, tmp1, tmp2, S1::T1, S2::T2) = GAPAData{T1,T2}(α12, tmp1, tmp2, S1, S2)

function init_algorithm!(alg::GAPA, model::AbstractFOSModel)
    S1, S2, n, status_generator = get_sets_and_status(alg, model)
    data = GAPAData(2.0, Array{Float64,1}(undef, n), Array{Float64,1}(undef, n),
            Ref(false), S1, S2)
    return data, status_generator
end

""" normedScalar(x1,x2,y1,y2)
Calculate |<x1-x2,y1-y2>|/(||x1-x2||*||y1-y2||)"""
function normedScalar(x1,x2,y1,y2)
    sum, norm1, norm2 = 0.0, 0.0, 0.0
    #@inbounds @simd
    for i in 1:length(x1)
        d1 = x1[i] - x2[i]
        d2 = y1[i] - y2[i]
        sum   += d1*d2
        norm1 += d1^2
        norm2 += d2^2
    end
    return abs(sum)/sqrt(norm1*norm2)
end

function S1constant!(y, alg::GAPA, data::GAPAData, mem, i)
    α12 = data.α12
    if !data.constinit.x
        y .= 0.0
        #TODO high accuracy?
        prox!(mem, S1, y)
        data.constinit.x = true
    end
    y .= α12.*mem# .+ (1-α1).*0.0
    return
end

function S1!(y, alg::GAPA, data::GAPAData, x, i, status::AbstractStatus, longstep=nothing)
    α12, S1 =  data.α12, data.S1
    prox!(y, S1, x)
    addprojeq(longstep, y, x)
    y .= α12.*y .+ (1-α12).*x
end

function S2!(y, alg::GAPA, data::GAPAData, x, i, status::AbstractStatus, longstep=nothing)
    α12, S2 =  data.α12, data.S2
    prox!(y, S2, x)
    checkstatus(status, y)
    addprojineq(longstep, y, x)
    y .= α12.*y .+ (1-α12).*x
end

function Base.step(alg::GAPA, data::GAPAData, x, i, status::AbstractStatus, longstep=nothing)
    α,β = alg.α,alg.β
    α12,tmp1,tmp2,S1,S2 = data.α12,data.tmp1,data.tmp2,data.S1,data.S2
    # Relaxed projection onto S1
    S1!(tmp1, alg, data, x,    i, status, longstep)
    # prox!(tmp1, S1, x)
    # addprojeq(longstep, tmp1, x)
    # tmp1 .= α12.*tmp1 .+ (1-α12).*x
    # Relaxed projection onto S2
    S2!(tmp2, alg, data, tmp1, i, status, longstep)
    # prox!(tmp2, S2, tmp1)
    # checkstatus(status, tmp2)
    # #println("α12: $α12")
    # addprojineq(longstep, tmp2, tmp1)
    # tmp2 .= α12.*tmp2 .+ (1-α12).*tmp1
    #Calculate θ_F estimate, normedScalar(x1,x2,y1,y2) = |<x1-x2,y1-y2>|/(||x1-x2||*||y1-y2||)
    scl = clamp(normedScalar(tmp2,tmp1,tmp1,x), 0.0, 1.0)
    s = sqrt(1-scl^2) #Efficient way to sin(acos(scl))
    # Update α1, α2
    αopt = 2/(1+s) # α1,α2 = 2/(1+sin(Θ_F))
    data.α12 = (1-β)*αopt + β*2.0
    #Set output
    x .= α.*tmp2 .+ (1-α).*x
    return
end

function getsol(::GAPA, data::GAPAData, x)
    tmp1,tmp2,S1,S2 = data.tmp1,data.tmp2,data.S1,data.S2
    prox!(tmp1, S1, x)
    prox!(tmp2, S2, tmp1)
    return tmp2
end

support_longstep(::GAPA) = true
projections_per_step(::GAPA) = (1,1)

support_linesearch(::GAPA) = Val{:Fast}
