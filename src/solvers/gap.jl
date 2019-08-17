
"""
GAP

"""
mutable struct GAP <: FOSAlgorithm
    α::Float64
    α1::Float64
    α2::Float64
    direct::Bool
    options
end
GAP(α=0.8, α1=1.8, α2=1.8; direct=false, kwargs...) = GAP(α, α1, α2, direct, kwargs)

struct GAPData{T1,T2} <: FOSSolverData
    tmp1::Array{Float64,1}
    tmp2::Array{Float64,1}
    constinit::Base.RefValue{Bool} # If constant term has been calculated (linesearch feature)
    S1::T1
    S2::T2
end

function init_algorithm!(alg::GAP, model::AbstractFOSModel)
    S1, S2, n, status_generator = get_sets_and_status(alg, model)
    data = GAPData(Array{Float64,1}(undef, n), Array{Float64,1}(undef, n),
            Ref(false), S1, S2)
    return data, status_generator
end

function S1constant!(y, alg::GAP, data::GAPData, mem, i)
    if !data.constinit.x
        y .= 0.0
        #TODO high accuracy?
        prox!(mem, S1, y)
        mem .= α1.*mem# .+ (1-α1).*0.0
        data.constinit.x = true
    end
    y .= mem
    return
end

function S1!(y, alg::GAP, data::GAPData, x, i, status::AbstractStatus, longstep=nothing)
    α1, S1 =  alg.α1, data.S1
    prox!(y, S1, x)
    addprojeq(longstep, y, x)
    y .= α1.*y .+ (1-α1).*x
end

function S2!(y, alg::GAP, data::GAPData, x, i, status::AbstractStatus, longstep=nothing)
    α2, S2 =  alg.α2, data.S2
    prox!(y, S2, x)
    checkstatus(status, y)
    addprojineq(longstep, y, x)
    y .= α2.*y .+ (1-α2).*x
end

function Base.step(alg::GAP, data::GAPData, x, i, status::AbstractStatus, longstep=nothing)
    α = alg.α
    tmp1, tmp2 = data.tmp1, data.tmp2
    # Relaxed projection onto S1
    S1!(tmp1, alg, data, x, i, status, longstep)
    # α1, S1 =  alg.α1, data.S1
    # prox!(tmp1, S1, x)
    # addprojeq(longstep, tmp1, x)
    # tmp1 .= α1.*tmp1 .+ (1-α1).*x
    #Relaxed projection onto S2
    S2!(tmp2, alg, data, tmp1, i, status, longstep)
    # α2, S2 =  alg.α2, data.S2
    # prox!(tmp2, S2, tmp1)
    # checkstatus(status, tmp2)
    # addprojineq(longstep, tmp2, tmp1)
    # tmp2 .= α2.*tmp2 .+ (1-α2).*tmp1
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

support_linesearch(::GAP) = Val{:Fast}

support_longstep(::GAP) = true
projections_per_step(::GAP) = (1,1)
