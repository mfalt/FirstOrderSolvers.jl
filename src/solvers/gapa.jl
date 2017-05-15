"""
GAP Adaptive

"""
type GAPA  <: FOSAlgorithm
    options
end
GAPA{T1,T2}(; kwargs...) = GAPA(kwargs)

type GAPAData{T1,T2} <: FOSSolverData
    α12::Float64
    tmp1::Array{Float64,1}
    tmp2::Array{Float64,1}
    S1::T1
    S2::T2
end
GAPAData{T1,T2}() = GAPAData{T1,T2}(α12, similar(x), similar(x), S1::T1, S2::T2)

function init_algorithm(::GAPA, model::FOSMathProgModel)
    hsde = HSDE(model)
    GAPAData(1.0, Array{Float64,1}(hsde.n), Array{Float64,1}(hsde.n), hsde.indAffine, hsde.indCones)
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

function Base.step(alg::GAPA, data::GAPAData, x, longstep=nothing)
    α12,tmp1,tmp2,S1,S2 = data.α12,data.tmp1,data.tmp2,data.S1,data.S2
    # Relaxed projection onto S1
    prox!(tmp1, S1, x)
    addprojeq(longstep, tmp1, x)
    tmp1 .= α12.*tmp1 .+ (1-α12).*x
    # Relaxed projection onto S2
    prox!(tmp2, S2, tmp1)
    addprojineq(longstep, tmp2, tmp1)
    tmp2 .= α12.*tmp2 .+ (1-α12).*tmp1
    #Calculate θ_F estimate, normedScalar(x1,x2,y1,y2) = |<x1-x2,y1-y2>|/(||x1-x2||*||y1-y2||)
    scl = clamp(normedScalar(tmp2,tmp1,tmp1,x), 0.0, 1.0)
    s = sqrt(1-scl^2) #Efficient way to sin(acos(scl))
    # Update α1, α2
    data.α12 = 2/(1+s) # α1,α2 = 2/(1+sin(Θ_F))
    #Set output
    x .= tmp2
    return
end

support_longstep(Alg::GAPA) = true
projections_per_step(Alg::GAPA) = (1,1)
