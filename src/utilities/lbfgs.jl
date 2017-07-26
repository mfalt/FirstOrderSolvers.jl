#Implementation based on Numerical Optimization, Jorge Nocedal
immutable LBFGSstate
    s::Array{Float64,2}
    y::Array{Float64,2}
    α::Array{Float64,1}
    ρ::Array{Float64,1}
    xold::Array{Float64,1}      #Temporary storage for xk
    gold::Array{Float64,1}      #Temporary storage for ∇f(xk)
    γ::Base.RefValue{Float64}   # H0 = γI = s'y/(y'y)*I if usegamma
    usegamma::Bool
    n::Int64                    #Total number possible vectors
    j::Base.RefValue{Int64}     #Location of next pair to be replaced
end

"""
LBFGSstate(m,n), state for using LBFGS
 m = length of each vector
 n = number of vectors to store
"""
LBFGSstate(m,n; usegamma=true) = LBFGSstate(zeros(m,n), zeros(m,n), zeros(n),
                                            zeros(n), zeros(m), zeros(m),
                                            Ref(1.0), usegamma,
                                            n, Ref(1))


 """
 initlbfgs!(x, g, state::LBFGSstate)

 Initialize LBFGSstate with new value x=x^{k+1} and g = ∇f(x^{k+1})
 """
 function initlbfgs!(x, g, state::LBFGSstate)
     copy!(state.xold, x)
     copy!(state.gold, g)
     state.j.x = state.j.x + 1
     return
 end

"""
set_x_∇f!(x, g, state::LBFGSstate)

Update LBFGSstate with new value x=x^{k+1} and g = ∇f(x^{k+1})
"""
function set_x_∇f!(x, g, state::LBFGSstate)
    @assert length(x) == length(state.xold)
    @assert length(g) == length(state.xold)

    j = state.j.x
    k = mod(j-1,state.n)+1
    state.s[:,k] .= x .- state.xold   # x_{k+1}) - x_{k}
    state.y[:,k] .= g .- state.gold   # ∇f(x_{k+1}) - ∇f(x_{k})
    copy!(state.xold, x)
    copy!(state.gold, g)
    tmp = vdot(state.s,state.y,k) #tmp = s[:,k]'y[:,k]'
    state.ρ[k] = 1/tmp
    state.j.x = j + 1

    if state.usegamma
        tmp = vdot(state.s,state.y,k)/vdot(state.y,state.y,k)

        if tmp < 1e6 && tmp > -1e6 #Some safety check
            state.γ.x = tmp
        end
        #println(state.γ.x)
    end
    return
end

"""
vdot(x::AbstractMatrix,y::AbstractVector,i) = dot(x[:,i],y)
vdot(x::AbstractMatrix,y::AbstractMatrix,i) = dot(x[:,i],y[:,i])
"""
function vdot(x::AbstractMatrix,y::AbstractVector,i)
    s = 0.0
    @inbounds @simd for j = 1:length(y)
        s += x[j,i]*y[j]
    end
    return s
end

function vdot(x::AbstractMatrix,y::AbstractMatrix,i)
    s = 0.0
    @inbounds @simd for j = 1:size(y,1)
        s += x[j,i]*y[j,i]
    end
    return s
end

"""
H∇f!(q, state::LBFGSstate)

Calculate q .= Hk*q given the LBFGSstate
"""
function H∇f!(q, state::LBFGSstate)
    s,y,α,ρ,n,j = state.s,state.y,state.α,state.ρ,state.n,state.j.x
    @assert length(q) == length(state.xold)
    #Backwards
    @inbounds for i = (n+j-1):-1:j
        k = mod(i-1,n)+1
        α[k] = ρ[k]*vdot(s,q,k)             # ρ[k]*s[:,k]'q
        LinAlg.axpy!(-α[k], view(y,:,k), q) # q .= q .- α[k].*y[:,k]
    end

    scale!(q, state.γ.x) # q = γ*q, i.e. H0 = γ*I
    #q .= q.*state.γ.x
    #Forwards
    @inbounds for i = j:(n+j-1)
        k = mod(i-1,n)+1
        β = ρ[k]*vdot(y,q,k)                 # β = ρ[k]*y[:,k]'q
        LinAlg.axpy!(α[k]-β, view(s,:,k), q) # q .= q .+ (α[k]-β)*s[:,k]
    end
    return #Result H*q is stored in q
end
