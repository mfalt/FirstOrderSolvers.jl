
"""
SuperMann

"""
type SuperMann <: FOSAlgorithm
    α::Float64
    c0::Float64
    c1::Float64
    q::Float64
    β::Float64
    σ::Float64
    λ::Float64
    direct::Bool
    nstore::Int64        #Number of vectors to store in LBFGS
    options
end

function SuperMann(c0=0., c1=0.99, q=0.99, β=0.5, σ=0.01, λ=1.0, nstore = 30; direct=false, kwargs...)
    α = 0.5 #TODO
    SuperMann(α,c0,c1,q,β,σ,λ,nstore,direct,kwargs)
end

immutable SuperMannData{T1,T2} <: FOSSolverData
    w::Array{Float64,1}
    d::Array{Float64,1}
    Rx::Array{Float64,1}
    Rw::Array{Float64,1}
    tmp1::Array{Float64,1}
    tmp2::Array{Float64,1}
    lbfgs::LBFGSstate
    S1::T1
    S2::T2
end

function init_algorithm!(alg::SuperMann, model::FOSMathProgModel)
    hsde, status_generator = HSDE(model, direct=alg.direct)
    data = SuperMannData(Array{Float64,1}(hsde.n), Array{Float64,1}(hsde.n),
                         Array{Float64,1}(hsde.n), Array{Float64,1}(hsde.n),
                         Array{Float64,1}(hsde.n), Array{Float64,1}(hsde.n),
                         LBFGSstate(hsde.n, alg.nstore),
                         hsde.indAffine, hsde.indCones)
    return data, status_generator
end

function S1!(y, alg::SuperMann, data::SuperMannData, x, i, status::AbstractStatus, longstep=nothing)
    α1 = 1.0
    prox!(y, data.S1, x)
    y .= α1.*y .+ (1-α1).*x
end

function S2!(y, alg::SuperMann, data::SuperMannData, x, i, status::AbstractStatus, longstep=nothing)
    α2 = 1.0
    prox!(y, data.S2, x)
    checkstatus(status, y)
    y .= α2.*y .+ (1-α2).*x
end

function R!(Rx, x, alg, data, i, status)
    S1!(tmp1, alg, data, x,    i, status, nothing)
    S2!(tmp2, alg, data, tmp1, i, status, nothing)
    #TODO Averaging?
    Rx .= x .- tmp2 # Rx = (I-T)x
    return
end


function Base.step(alg::SuperMann, data::SuperMannData, x, i, status::AbstractStatus, longstep=nothing)
    α, c0, c1, q, β, σ, λ = alg.α, alg.c0, alg.c1, alg.q, alg.β, alg.σ, alg.λ
    w, d, Rx, Rw, lbfgs = data.w, data.d, data.x, data.w, data.lbfgs
    tmp1, tmp2 = data.tmp1, data.tmp2
    #ηk    = data.η.x     #Reference
    rsafe = data.rsafe.x  #Reference

    R!(Rx, x, alg, data, i, status)     # Rx = R(x), and do status check

    isdone(status) && return            # Step 1, return on convergence

    set_x_∇f!(x, Rx, lbfgs)             # Update the lbfgs to remember these values
    normRx = norm(Rx)

    #Calculate newton direction d:
    copy!(d,Rx)    # d .= Rx
    H∇f!(d, lbfgs) # d .= H*d             # Step 2

    if normRx <= c0*data.η.x            # Step 3
        data.η.x = normRx               # Step 3
        x .= x .+ d                     # Blind update
    else
        # nk = nk                       # Step 4
        τk = 1.0                        # Step 4
        while true
            w .= x .+ τk.*d             # Step 5
            R!(Rw, w, alg, data, i, NoStatus(:no)) #Rw = R(w)
            normRw = norm(Rw)
            if normRx <= rsafe && normRw <= c1*normRx
                x .= w                  # Step 5a
                rsafe = normRx + q      # Step 5a
                break                   # Go to step 6
            else
                ρk = normRw^2 - 2α*sum(Rw.*(w.-x)) #dot(Rw, w-x)
                if ρk >= σ*normRw*normRx
                    scal = λ*ρk/normRw^2
                    x .= x .- scal.*Rw  # Step 5b
                else
                    τk = β*τk           # Step 5b otherwise
                end                     # Return Step 5
            end
        end
    end
    data.rsafe.x = rsafe
    return
end

function getsol(::SuperMann, data::SuperMannData, x)
    tmp1,tmp2,S1,S2 = data.tmp1,data.tmp2,data.S1,data.S2
    prox!(tmp1, S1, x)
    prox!(tmp2, S2, tmp1)
    return tmp2
end
