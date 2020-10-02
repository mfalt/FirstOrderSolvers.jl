export AndersonWrapper

mutable struct AndersonWrapper{T<:FOSAlgorithm} <: FOSAlgorithm
    lsinterval::Int64
    alg::T
    options
end

mutable struct AndersonWrapperData{T<:FOSSolverData} <: FOSSolverData
    lsinterval::Int64
    tmp1::Array{Float64,1}
    tmp2::Array{Float64,1}
    tmp3::Array{Float64,1}
    mem::Array{Float64,1}
    res::Array{Float64,1}
    algdata::T
end

function AndersonWrapper(alg::T; lsinterval=100, kwargs...) where T
    if support_linesearch(alg) == Val{:False}
        @error "Algorithm $T does not support line search"
    end
    AndersonWrapper{T}(lsinterval, alg, merge(alg.options, kwargs))
end

function init_algorithm!(ls::AndersonWrapper, model::AbstractFOSModel)
    alg, lsinterval = ls.alg, ls.lsinterval
    algdata, status =  init_algorithm!(alg, model)
    x = getinitialvalue(model, alg, algdata)
    data = AndersonWrapperData(lsinterval, similar(x), similar(x), similar(x), similar(x), similar(x), algdata)
    return data, status
end

#getmn(data::AndersonWrapperData) = getmn(data.algdata)

function ϕ(θb, η)
    if η >= θb
        return one(θb)
    else
        sη = η<0 ? -1 : 1 # Sign with 0 -> 1
        return (1-sn*θb)/(1-η)
    end
end

function Base.step(wrap::AndersonWrapper, adata::AndersonWrapperData, xk, i, status::AbstractStatus, longstep=nothing)

    gxk, gxk1 # g(xᵏ) :=xᵏ-f(xᵏ), g(xᵏ⁻¹)
    gx̃k # g(x̃ᵏ)
    yk1 # yₖ₋₁
    s[:,i] # Array with sₖ, where i=mk -> k-1
    ŝ[:,i] # # Array with ŝₖ, where i=mk -> k-1
    x̃k # x̃ᵏ
    xk1 # xᵏ⁻¹
    Hỹ[:,i] # Array with Hₖ₋ⱼỹₖ₋ⱼ, where i=mk -> k-1
    Hs[:,i] # Array with Hₖ₋ⱼ'ŝₖ₋ⱼ/(ŝₖ₋ⱼ'Hyₖ₋ⱼ), where i=mk -> k-1

    # Constants:
    θb # θ̄  Regularization
    τ  # Regularization
    # α ? 
    D # Safe-guard
    ϵ # Safe-guard
    Ū # Inital residual
    Na # N_{AA}
    m # Max memory

    adata.mk = adata.mk + 1 # Line 4
    mk = adata.mk # For local use

    s[:,mk] .= x̃k .- xk1        # Line 5
    yk1 .= gx̃k .- gxk1      # Line 5
    # ŝk1 .= sk1 .- sum(ŝj -> (ŝj'sk1)/(ŝj'ŝj.*ŝj), ŝ[(k-mk):k-2])
    # Line 6
    ŝ[:,mk] .= s[:,mk] .+ sum(j -> ((ŝ[j]'s[:,mk])/(ŝ[:,j]'ŝ[:,j])).*ŝ[:,j] , 1:(mk-1))
    # Line 7-8
    if mk == m+1 || norm(ŝk1) < τ*norm(sk1)
        adata.mk = 1
        mk = 1
        ŝ[:,mk] .= s[:,mk]
        # Hk1 = I
    end
    # Line  9-10
    # TODO use Hy = Hk1*yk1 = Hk1*(gx̃-gxk1)
    γk1 = ŝ[:,mk]'Hy[:,mk]/norm(ŝ[:,mk])^2
    θk1 = ϕ(θb, γk1)
    ỹk1 .= θk1.*yk1 .- (1-θk1).*gxk1

    # TODO Compute Hỹ
    # TODO Compute Hs

    # Line 11
    # Hk1 .= Hk1 + ... # Next Hk1
    # x̃k .= xk .- Hk1*gxk # Next x̃k
    dk .= gxk + sum(j -> (s[:,j] .- Hỹ[:,j]).* (Hs[:,j]'gxk), 1:mk)
    x̃k .= xk .- dk  # Next index x̃

    # Line 12-14
    if norm(gxk) <= D*U/(Na +1)^(1+ϵ)
        xk1 .= xk # Save for next index
        xk .= x̃k # Next index x and x̃

        # We need to compute f(xk) now, use xk as temp variable for f(xk)
        step(wrap.alg, adata.algdata, xk, i, status, longstep) # We only use this for gxk, gx̃k
        gxk1 .= gx̃k # Save for next index

        # Compute for next index
        gxk .= x̃k .- xk  # this is same as x̃k - f(x̃k)
        gx̃k .= gxk
       
        adata.Na += 1
    else
        xk1 .= xk   # Save for next index
        step(wrap.alg, adata.algdata, xk, i, status, longstep) # Next index xk
        gxk1 .= gxk # Save for next index
        gxk .= xk1 .- xk # Next index
    end
end

function normdiff(x::T,y::T) where T<:StridedVector
    nx, ny = length(x), length(y)
    s = 0.0
    @assert nx == ny
    for j = 1:nx
        s += abs2(x[j]-y[j])
    end
    return sqrt(s)
end

function get_normres!(wrap::AndersonWrapper, adata::AndersonWrapperData, x, i, status::AbstractStatus)
    adata.tmp2 .= x
    step(wrap.alg, adata.algdata, x, i, status, nothing)
    #If time to do adata
    normres = normdiff(x,adata.tmp2)
end

function getsol(alg::AndersonWrapper, data::AndersonWrapperData, x)
    getsol(alg.alg, data.algdata, x)
end
