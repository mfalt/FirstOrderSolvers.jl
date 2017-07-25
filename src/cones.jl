import ProximalOperators: ProximableFunction, prox!, IndPSD

# TODO zero cone
const conemap = Dict{Symbol, ProximableFunction}(
    :Free => IndFree(),
    :Zero => IndZero(),
    :NonNeg => IndNonnegative(),
    :NonPos => IndNonpositive(),
    :SOC => IndSOC(),
    :SOCRotated => IndRotatedSOC(),
    :SDP => IndPSD(scaling=true),
    :ExpPrimal => IndExpPrimal(),
    :ExpDual => IndExpDual()
)
const badcones = ProximableFunction[]


#Wrapper for dual prox to avoid double duals
function proxDual!(y::AbstractArray, C::ProximableFunction, x::AbstractArray)
    prox!(y, C, -x)
    @simd for i = 1:length(x)
        y[i] = x[i] + y[i]
    end
end

## Some better dual projections
const myIndFree = IndFree()
proxDual!(y::AbstractArray, C::IndZero, x::AbstractArray) = prox!(y, myIndFree, x)
const myIndZero = IndPoint()
proxDual!(y::AbstractArray, C::IndFree, x::AbstractArray)          = prox!(y, myIndZero, x)
proxDual!(y::AbstractArray, C::IndNonnegative, x::AbstractArray)   = prox!(y, C, x)
proxDual!(y::AbstractArray, C::IndNonpositive, x::AbstractArray)   = prox!(y, C, x)
# TODO figure out if self dual PSD
#proxDual!(y::AbstractArray, C::IndPSD, x::AbstractArray)           = prox!(y, C, x)

# Unroll the loop over the different types of functions to prox on, same as for prox! in ProximalOperators
# Ignores function value
@generated function proxDual!{T, A, B, N}(y::AbstractArray{T}, f::SlicedSeparableSum{A, B, N}, x::AbstractArray{T})
  ex = :()
  for i = 1:N # For each function type
    ex = quote $ex;
      for k in eachindex(f.fs[$i]) # For each function of that type
        proxDual!(view(y, f.idxs[$i][k]...), f.fs[$i][k], view(x,f.idxs[$i][k]...))
      end
    end
  end
  ex = :($ex; return)
end

""" Indicator of K2×K1*×R+ × K2*×K1×R+ ∈ (R^n,R^m,R)^2"""
type DualConeProduct{T1<:SlicedSeparableSum,T2<:SlicedSeparableSum} <: ProximableFunction
    K1::T1
    K2::T2
    m::Int64
    n::Int64
end

"""
Finds the largest index in a SlicedSeparableSum
"""
largestindex(e) = error("Something unexpected happend with the product cones with index of type $(typeof(e))")
largestindex(K::SlicedSeparableSum) = largestindex(K.idxs)                      # Recursive call on indices
largestindex(l::Union{Tuple, AbstractArray})     = maximum(largestindex.(l))    # Recursive call, apply on each element
largestindex(l::Union{AbstractUnitRange,Number}) = maximum(l)                   # Possible elements

DualConeProduct{T1,T2}(K1::T1,K2::T2) = DualConeProduct{T1,T2}(K1, K2, largestindex(K1), largestindex(K2))

function prox!(y::AbstractArray, K::DualConeProduct, x::AbstractArray)
    m, n = K.m, K.n
    K1, K2 = K.K1, K.K2
    nu = n+m+1
    xx = view(x, 1:n)
    xy = view(x, (n+1):(n+m))
    xr = view(x, (nu+1):(nu+n))
    xs = view(x, (nu+n+1):(nu+n+m))

    yx = view(y, 1:n)
    yy = view(y, (n+1):(n+m))
    yr = view(y, (nu+1):(nu+n))
    ys = view(y, (nu+n+1):(nu+n+m))

    prox!(    yx, K2, xx)
    proxDual!(yy, K1, xy)
    y[nu] = max(x[nu], 0)
    proxDual!(yr, K2, xr)
    prox!(    ys, K1, xs)
    y[2nu] = max(x[2nu], 0)
end
