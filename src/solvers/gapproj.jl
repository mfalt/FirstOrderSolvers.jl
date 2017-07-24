
"""
GAPP projected GAP
Defaults to direct solve of the linear system
"""
type GAPP <: FOSAlgorithm
    α::Float64
    α1::Float64
    α2::Float64
    iproj::Int64
    direct::Bool
    options
end
GAPP(α=0.8, α1=1.8, α2=1.8; direct=true, iproj=100, kwargs...) = GAPP(α, α1, α2, iproj, direct, kwargs)

immutable GAPPData{T1,T2} <: FOSSolverData
    tmp1::Array{Float64,1}
    tmp2::Array{Float64,1}
    S1::T1
    S2::T2
end

function init_algorithm!(alg::GAPP, model::FOSMathProgModel)
    hsde, status_generator = HSDE(model, direct=alg.direct)
    data = GAPPData(Array{Float64,1}(hsde.n), Array{Float64,1}(hsde.n), hsde.indAffine, hsde.indCones)
    return data, status_generator
end

function Base.step(alg::GAPP, data::GAPPData, x, i, status::AbstractStatus, longstep=nothing)
    α,α1,α2 = alg.α, alg.α1, alg.α2
    tmp1,tmp2,S1,S2 = data.tmp1, data.tmp2, data.S1, data.S2
    # Relaxed projection onto S1
    prox!(tmp1, S1, x)
    if i % alg.iproj == 0
        tmp3 = similar(tmp1)
        tmp4 = similar(tmp1)
        res = similar(tmp1)
        # Projection onto S2
        prox!(tmp2, S2, tmp1)
        prox!(res, S1, tmp2)
        res .= res - tmp1
        #res = Ps1*(Ps2*Ps1-I)x

        normbest = Inf
        αbest = -1.0
        for k = 0:20
            αtest = 2.0.^k
            tmp3 .= tmp1 .+ αtest.*res
            prox!(tmp4, S2, tmp3)
            normtest = norm(tmp4 - tmp3)
            println("normtest: $normtest")
            if normtest < normbest
                αbest = αtest
                normbest = normtest
            end
        end
        println("αbest: $αbest")
        tmp1 .= tmp1 .+ αbest.*res
        prox!(tmp2, S2, tmp1)
        checkstatus(status, tmp2)
        tmp2 .= α2.*tmp2 .+ (1-α2).*tmp1
        x .= tmp2
    else
        tmp1 .= α1.*tmp1 .+ (1-α1).*x
        # Relaxed projection onto S2
        prox!(tmp2, S2, tmp1)
        checkstatus(status, tmp2)
        tmp2 .= α2.*tmp2 .+ (1-α2).*tmp1
        # Relaxation
        x .= α.*tmp2 .+ (1-α).*x
    end

    return
end

function getsol(::GAPP, data::GAPPData, x)
    tmp1,tmp2,S1,S2 = data.tmp1,data.tmp2,data.S1,data.S2
    prox!(tmp1, S1, x)
    prox!(tmp2, S2, tmp1)
    return tmp2
end

support_longstep(alg::GAPP) = false
projections_per_step(alg::GAPP) = (0,0)
