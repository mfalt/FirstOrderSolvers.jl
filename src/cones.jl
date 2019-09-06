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

# type DualCone{T<:ProximableFunction} <: ProximableFunction
#     C::T
# end

# dual{T<:ProximableFunction}(C::T) = DualCone{T}(C)
# dual(C::DualCone) = C.C

# function prox!(cone::DualCone, x, y)
#     prox!(cone.C, -x, y)
#     for i = 1:length(x)
#         y[i] = x[i] + y[i]
#     end
# end

mutable struct ConeProduct{N,T} <: ProximableFunction
    ranges::NTuple{N, UnitRange{Int64}}
    cones::T
end
#
ConeProduct() = ConeProduct{0,Any}((),())

# function ConeProduct(ranges::NTuple{N, UnitRange{Int64}}, cones::T) where {N,T}
#     return ConeProduct{N,T}(ranges, cones)
# end

toRanges(rangesIn::NTuple{N, UnitRange{Int64}}) where N = rangesIn

function toRanges(rangesIn::NTuple{N,Array{T,1}}) where {N, T<:Integer}
    ranges = Array{UnitRange{Int64},1}(undef, N)
    for j in 1:N
        range = rangesIn[j]
        for (i,el) in enumerate(range[1]:range[end])
            if el != range[i]
                @error "Invalid range in input"
            end
        end
        ranges[j] = range[1]:range[end]
    end
    return tuple(ranges...)
end

#Should accept tro tuples with UnitRange{Int64} and ProximableFunction
#Cant enforce types or ambigous with default constructor?
function ConeProduct(rangesIn, cones)
    ranges = toRanges(rangesIn)
    N = length(cones)
    @assert typeof(ranges) == NTuple{N, UnitRange{Int64}}
    @assert length(ranges) == N
    @assert typeof(cones) <: Tuple
    #Verify that ranges covers the whole range
    prevRangeEnd = 0
    for i = 1:N
      @assert ranges[i][1] == prevRangeEnd + 1
      prevRangeEnd = ranges[i][end]
      @assert typeof(cones[i]) <: ProximableFunction
    end
    T = typeof(cones)
    #println("Dumping cones input")
    #dump(cones)
    ConeProduct(ranges, cones)
end

#Wrapper for dual prox to avoid double duals
function proxDual!(y::AbstractArray, C::ProximableFunction, x::AbstractArray)
    prox!(y, C, -x)
    @simd for i = 1:length(x)
        y[i] = x[i] + y[i]
    end
end
# proxDual!(y, C::DualCone, x) = prox!(y, C.C, x)


function prox!(y::AbstractArray, C::ConeProduct{N,T}, x::AbstractArray) where {N,T}
    #TODO Paralell implementation
    for i = 1:N
        prox!(view(y, C.ranges[i]), C.cones[i], view(x, C.ranges[i]))
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

function proxDual!(y::AbstractArray, C::ConeProduct{N,T}, x::AbstractArray) where {N,T}
    #TODO Paralell implementation
    for i = 1:N
        proxDual!(view(y, C.ranges[i]), C.cones[i], view(x, C.ranges[i]))
    end
end

""" Indicator of K2×K1*×R+ × K2*×K1×R+ ∈ (R^n,R^m,R)^2"""
mutable struct DualConeProduct{T1<:ConeProduct,T2<:ConeProduct} <: ProximableFunction
    K1::T1
    K2::T2
    m::Int64
    n::Int64
end

DualConeProduct(K1::T1,K2::T2) where {T1,T2} = DualConeProduct{T1,T2}(K1, K2, K1.ranges[end][end], K2.ranges[end][end])
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
