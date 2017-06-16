
"""
GAP

"""
type GAP <: FOSAlgorithm
    α::Float64
    α1::Float64
    α2::Float64
    options
end
GAP(α=0.8, α1=1.8, α2=1.8; kwargs...) = GAP(α, α1, α2, kwargs)

immutable GAPData{T1,T2} <: FOSSolverData
    tmp1::Array{Float64,1}
    tmp2::Array{Float64,1}
    S1::T1
    S2::T2
end

function init_algorithm!(::GAP, model::FOSMathProgModel)
    hsde = HSDE(model)
    GAPData(Array{Float64,1}(hsde.n), Array{Float64,1}(hsde.n), hsde.indAffine, hsde.indCones)
end

function Base.step(alg::GAP, data::GAPData, x, i, status::AbstractStatus, longstep=nothing)
    α,α1,α2 = alg.α, alg.α1, alg.α2
    tmp1,tmp2,S1,S2 = data.tmp1, data.tmp2, data.S1, data.S2
    # Relaxed projection onto S1
    prox!(tmp1, S1, x)
    addprojeq(longstep, tmp1, x)
    tmp1 .= α1.*tmp1 .+ (1-α1).*x
    # Relaxed projection onto S2
    prox!(tmp2, S2, tmp1)
    checkstatus(status, tmp2)
    addprojineq(longstep, tmp2, tmp1)
    tmp2 .= α2.*tmp2 .+ (1-α2).*tmp1
    # Relaxation
    x .= α.*tmp2 .+ (1-α).*x
    return
end

function getsol(::GAP, data::GAPData, x)
    tmp1,tmp2,S1,S2 = data.tmp1,data.tmp2,data.S1,data.S2
    prox!(tmp1, S1, x)
    prox!(tmp2, S2, tmp1)
    return tmp2
end

support_longstep(::GAP) = true
projections_per_step(::GAP) = (1,1)
