export AndersonWrapper

mutable struct AndersonWrapper{T<:FOSAlgorithm} <: AbstractWrapper
    alg::T
    α::Float64
    m::Int64 # Max memory
    θb::Float64 # θ̄  Regularization
    τ::Float64  # Regularization
    D::Float64 # Safe-guard
    ϵ::Float64 # Safe-guard
    options
end

mutable struct AndersonWrapperData{T<:FOSSolverData} <: AbstractWrapperData
    algdata::T

    # Algorithm memory
    xk::Array{Float64,1}
    gxk::Array{Float64,1}
    gxk1::Array{Float64,1}
    gx̃k::Array{Float64,1}
    yk1::Array{Float64,1}
    ỹk1::Array{Float64,1}
    s::Array{Float64,2}
    ŝ::Array{Float64,2}
    ŝk1::Array{Float64,1}
    x̃k::Array{Float64,1}
    xk1::Array{Float64,1}
    Hỹ::Array{Float64,2}
    Hs::Array{Float64,2}
    Hy::Array{Float64,1}
    x0::Array{Float64,1}
    mk::Int # Max memory

    # Constants:
    # α ? 
    Ū::Float64 # Inital residual
    Na::Int # N_{AA}
end

function AndersonWrapper(alg::T; α=0.9, m=10, θ=0.01, τ=0.001, D = 1e6, ϵ=1e-6, kwargs...) where T
    # if support_linesearch(alg) == Val{:False}
    #     @error "Algorithm $T does not support line search"
    # end
    AndersonWrapper{T}(alg, α, m, θ, τ, D, ϵ, merge(alg.options, kwargs))
end

function init_algorithm!(aa::AndersonWrapper, model::AbstractFOSModel)
    algdata, status =  init_algorithm!(aa.alg, model)
    x = getinitialvalue(model, aa.alg, algdata)
    nx, m = length(x), aa.m
    aadata = AndersonWrapperData(algdata, similar(x), similar(x), similar(x), similar(x), similar(x), similar(x),
        similar(x, nx, m+1), similar(x, nx, m+1), similar(x), similar(x), similar(x), similar(x, nx, m), similar(x, nx, m), similar(x),
        copy(x),
        0, 0.0, 0)
    return aadata, status
end

#getmn(data::AndersonWrapperData) = getmn(data.algdata)

function ϕ(θb, η)
    if η >= θb
        return one(θb)
    else
        sη = η<0 ? -1 : 1 # Sign with 0 -> 1
        return (1-sη*θb)/(1-η)
    end
end

""" H!(dk, gxk, s, Hỹ, Hs, mk)
    Compute dk = H*gxk
""" 
function H!(dk, gxk, s, Hỹ, Hs, mk)
    #dk .= gxk .+ sum(j -> (s[:,j] .- Hỹ[:,j]) .* (Hs[:,j]'gxk), 1:mk)
    dk .= gxk .+ (s[:,1:mk] - Hỹ[:,1:mk])*(Hs[:,1:mk]'gxk)
    return
end

function Base.step(wrap::AndersonWrapper, adata::AndersonWrapperData, xk, i, status::AbstractStatus, longstep=nothing)
    gxk = adata.gxk # g(xᵏ) :=xᵏ-f(xᵏ)
    gxk1 = adata.gxk1 # g(xᵏ⁻¹)
    gx̃k = adata.gx̃k# g(x̃ᵏ)
    yk1 = adata.yk1 # yₖ₋₁
    ỹk1 = adata.ỹk1 # ỹₖ₋₁
    s = adata.s # Array with sₖ, where i=mk -> k-1
    ŝ = adata.ŝ # # Array with ŝₖ, where i=mk -> k-1, normalized
    ŝk1 = adata.ŝk1 # ŝₖ₋₁, not normalized
    x̃k = adata.x̃k # x̃ᵏ
    xk1 = adata.xk1 # xᵏ⁻¹
    Hỹ = adata.Hỹ # Array with Hₖ₋ⱼỹₖ₋ⱼ, where i=mk -> k-1
    Hs = adata.Hs # Array with Hₖ₋ⱼ'ŝₖ₋ⱼ/(ŝₖ₋ⱼ'Hyₖ₋ⱼ), where i=mk -> k-1
    Hy = adata.Hy # Hₖ₋ⱼyₖ₋ⱼ

    Ū = adata.Ū # Inital residual
    Na = adata.Na # N_{AA}

    # Constant parameters:
    α = wrap.α # θ̄  Regularization
    θb = wrap.θb # θ̄  Regularization
    τ = wrap.τ  # Regularization
    # α ? 
    D = wrap.D # Safe-guard
    ϵ = wrap.ϵ # Safe-guard
    m = wrap.m # Max memory


    adata.mk = adata.mk + 1 # Line 4
    mk = adata.mk # For local use

    if i == 1
        x0 = adata.x0
        x1 = copy(x0)
        step(wrap.alg, adata.algdata, x1, 0, NoStatus(:Continue), nothing)
        adata.xk1 .= x0
        adata.gxk1 .= x0 .- x1
        adata.Ū = norm(adata.gxk1)
    
        adata.xk .= α.* x1 .+ (1-α) .* x0
        adata.x̃k .= adata.xk
        fx1 = copy(adata.xk)
        step(wrap.alg, adata.algdata, fx1, 0, NoStatus(:Continue), nothing)
        adata.gxk .= x1 .- fx1
        adata.gx̃k .= adata.gxk
    end
    # println("Iteration: ",i)
    #println("xk: ", xk)
    s[:,mk] .= x̃k .- xk1        # Line 5
    yk1 .= gx̃k .- gxk1      # Line 5
    # ŝk1 .= sk1 .- sum(ŝj -> (ŝj'sk1)/(ŝj'ŝj.*ŝj), ŝ[(k-mk):k-2])
    # Line 6
    #ŝ[:,mk] .= s[:,mk] .+ sum(j -> ((ŝ[:,j]'s[:,mk])/(ŝ[:,j]'ŝ[:,j])).*ŝ[:,j] , 1:(mk-1))
    ŝk1 .= s[:,mk] .+ ŝ[:,1:(mk-1)]*(ŝ[:,1:(mk-1)]'s[:,mk]) # The c

    # Line 7-8
    if mk == m+1 || norm(ŝk1) < τ*norm(s[:,mk])
        # println("Resetting mk")
        adata.mk = 1
        mk = 1
        ŝk1 .= s[:,mk]
        # Hk1 = I
    end

    # Normalizing
    ŝ[:,mk] .= ŝk1./(ŝk1'ŝk1)
    
    # Line  9-10
    # TODO use Hy = Hk1*yk1 = Hk1*(gx̃-gxk1)
    H!(Hy, yk1, s, Hỹ, Hs, mk-1)
    γk1 = ŝ[:,mk]'Hy
    θk1 = ϕ(θb, γk1)
    # println("θk1: ",θk1)
    ỹk1 .= θk1.*yk1 .- (1-θk1).*gxk1

    
    # Compute Hỹ (TODO Fix indices)
    H!(view(Hỹ,:,mk), ỹk1, s, Hỹ, Hs, mk-1)
    # Compute Hs (TODO Fix indices)
    H!(view(Hs,:,mk), view(ŝ,:,mk), s, Hỹ, Hs, mk-1)

    # println("Hỹmk:", sum(Hỹ[:,mk]))
    # println("Hsk:", sum(Hs[:,mk]))
    # println("smk:", sum(s[:,mk]))
    # println("Hy : ", sum(Hy))
    # Line 11
    # Hk1 .= Hk1 + ... # Next Hk1
    # x̃k .= xk .- Hk1*gxk # Next x̃k
    dk = ỹk1 # Just renaming, we dont need ỹk1 anymore
    H!(dk, gxk, s, Hỹ, Hs, mk) # dk .= gxk + sum(j -> (s[:,j] .- Hỹ[:,j]).* (Hs[:,j]'gxk), 1:mk)
    x̃k .= xk .- dk  # Next index x̃

    # Compute f(x̃k) (next index)
    gx̃k .= x̃k
    step(wrap.alg, adata.algdata, gx̃k, i, status, longstep) # We only use this for gxk, gx̃k
    gx̃k .= x̃k .- gx̃k # g(x) = x - f(x)

    # Line 12-14
    # println("gk ", norm(gxk), "guard: ",D*Ū/(Na +1)^(1+ϵ))
    if norm(gxk) <= D*Ū/(Na +1)^(1+ϵ) #|| i == 1
        # println("Acc at :", i)
        xk1 .= xk # Save for next index
        gxk1 .= gxk # Save for next index

        xk .= x̃k # Next index x and x̃
        gxk .= gx̃k # Next index gx
       
        adata.Na += 1
    else
        # TODO α relaxation here
        xk1 .= xk   # Save for next index
        gxk1 .= gxk # Save for next index

        # Compute f(xk)
        step(wrap.alg, adata.algdata, xk, i, NoStatus(:Continue), longstep) # Next index xk now
        gxk .= xk1 .- xk # Next index
        xk .= α.* xk1 .+ (1-α) .* xk # Relaxation
    end
    # println("xk:", xk)
end
