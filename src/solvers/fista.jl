
"""
FISTA

"""
mutable struct FISTA <: FOSAlgorithm
    α::Float64
    direct::Bool
    options
end
FISTA(α=1.0; direct=false, kwargs...) = FISTA(α, direct, kwargs)

struct FISTAData{T1,T2} <: FOSSolverData
    t::Base.RefValue{Float64}
    y::Array{Float64,1}
    xold::Array{Float64,1}
    tmp1::Array{Float64,1}
    S1::T1
    S2::T2
end

function init_algorithm!(alg::FISTA, model::AbstractFOSModel)
    S1, S2, n, status_generator = get_sets_and_status(alg, model)
    data = FISTAData(Ref(1.0), zeros(n), zeros(n), Array{Float64,1}(undef, n), S1, S2)
    return data, status_generator
end

function Base.step(alg::FISTA, data::FISTAData, x, i, status::AbstractStatus, longstep=nothing)
    α,y,xold,tmp1,S1,S2 = alg.α,data.y,data.xold,data.tmp1,data.S1,data.S2
    t = data.t #Reference
    if i == 1 #TODO this is ugly initializing hack
        y .= x
    end
    # Relaxed projection onto S1
    prox!(tmp1, S1, y)
    addprojeq(longstep, tmp1, y)
    tmp1 .= α.*tmp1 .+ (1-α).*y
    # Relaxed projection onto S2
    xold .= x
    prox!(x, S2, tmp1)
    checkstatus(status, x)
    addprojineq(longstep, x, tmp1)

    told = t.x
    t.x = (1+sqrt(1+4*t.x^2))/2
    y .= x .+ (told-1)./(t.x).*(x .- xold)
    return
end

function getsol(::FISTA, data::FISTAData, x)
    tmp1,S1,S2 = data.tmp1,data.S1,data.S2
    tmp2 = similar(x)
    prox!(tmp1, S1, x)
    prox!(tmp2, S2, tmp1)
    return tmp2
end

support_longstep(::FISTA) = true
projections_per_step(::FISTA) = (1,1)
