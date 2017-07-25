
immutable HSDEMatrixQ{T,M<:AbstractMatrix{T},V<:AbstractVector{T}} <: AbstractMatrix{T}
    A::M
    b::V
    c::V
    am::Int64
    an::Int64
end

"""
HSDEMatrixQ(A::M, b::V, c::V) where {T,M<:AbstractMatrix{T},V<:AbstractVector{T}}
"""
function HSDEMatrixQ(A::M, b::V, c::V) where {T,M<:AbstractMatrix{T},V<:AbstractVector{T}}
    am, an = size(A)
    @assert size(b,1) == am
    @assert size(c,1) == an
    HSDEMatrixQ{T,M,V}(A, b, c, am, an)
end

Base.size(Q::HSDEMatrixQ) = (Q.am+Q.an+1, Q.am+Q.an+1)

function Base.showarray(io::IO, Q::T, val::Bool) where {T<:HSDEMatrixQ}
    m,n = size(Q)
    println(io, "$(m)x$(m) $T:")
    println(io, "[0   A'  c;\n -A  0   b;\n -c' -b' 0]")
    println(io,"where A:")
    Base.showarray(io, Q.A, false)
    println(io,"\nb:")
    Base.showarray(io, Q.b, true)
    println(io,"\nc:")
    Base.showarray(io, Q.c, true)
end

# Q of size
#(n,n) (n,m) (n,1)
#(m,m) (m,m) (m,1)
#(1,n) (1,m) (1,1)
function LinAlg.A_mul_B!{T,M<:AbstractArray{T,2},V<:AbstractArray{T,1}}(Y::AbstractArray{T,1}, Q::HSDEMatrixQ{T,M,V}, B::AbstractArray{T,1})
    @assert size(Q) == (size(Y,1), size(B,1))
    y1 = view(Y, 1:Q.an)
    y2 = view(Y, (Q.an+1):(Q.an+Q.am))
    b1 = view(B, 1:Q.an)
    b2 = view(B, (Q.an+1):(Q.an+Q.am))
    b3 = B[Q.am+Q.an+1]
    A = Q.A
    b = Q.b
    c = Q.c
    At_mul_B!(y1, A, b2)
    A_mul_B!( y2, A, b1)
    #TODO maybe switch back for readability?
    @blas! y1 += b3*c # y1 .=  y1 .+ c.*b3
    @blas! y2 -= b3*b # y2 .= b.*b3 .- y2
    scale!(y2, -1)
    Y[Q.am+Q.an+1] = - dot(c, b1) - dot(b,b2)
    return Y
end

function LinAlg.At_mul_B!{T,M<:AbstractArray{T,2},V<:AbstractArray{T,1}}(Y::AbstractArray{T,1}, Q::HSDEMatrixQ{T,M,V}, B::AbstractArray{T,1})
    A_mul_B!(Y,Q,B)
    scale!(Y, -1)
    return Y
end


immutable HSDEMatrix{T,QType<:AbstractMatrix{T}} <: AbstractMatrix{T}
    Q::QType
    cgdata::CGdata{T}
end

"""
HSDEMatrix(Q::HSDEMatrixQ)
"""
HSDEMatrix{T, QType<:AbstractMatrix{T}}(Q::QType) = HSDEMatrix{T,QType}(Q, CGdata(2*size(Q)[1]))

"""
HSDEMatrix(A::M, b::V, c::V) where {T,M<:AbstractMatrix{T},V<:AbstractVector{T}}
"""
function HSDEMatrix(A::M, b::V, c::V) where {T,M<:AbstractMatrix{T},V<:AbstractVector{T}}
    Q = HSDEMatrixQ(A, b, c)
    HSDEMatrix(Q, CGdata(2*size(M.Q)[1]))
end

function Base.size(M::HSDEMatrix)
    m, n = size(M.Q)
    return (2m,2m)
end

function Base.showarray(io::IO, M::T, val::Bool) where {T<:HSDEMatrix}
    m,n = size(M)
    println(io, "$(m)x$(m) $T:")
    println(io, "[I  Q';\n Q -I]")
    println(io,"where Q:")
    Base.showarray(io, M.Q, false)
end

"""
    solve argmin_y{||x-y||₂}, s.t. Q*u==v, where [u;v] == x
"""
function prox!(y::AbstractVector, A::HSDEMatrix, x::AbstractVector)
    tol = size(A,2)*eps()
    max_iters = 1000
    cgdata = A.cgdata
    if cgdata.firstrun.x #Pointer, has this been initialized
        cgdata.xinit .= x #Works as guess since Q square
        cgdata.firstrun.x = false
    end
    #Since y is the initial guess, use previous solution
    y .= cgdata.xinit
    #solve [I Q';Q -I][u^(k+1);μ] = [u^k;v^k]
    conjugategradient!(y, A, x, cgdata.r, cgdata.p, cgdata.z,
                       tol = tol, max_iters = max_iters)
    #Save initial guess for next prox!
    cgdata.xinit .= y
    m, n = size(A.Q)
    #Let v=Q*u
    u = view(y,1:m)
    v = view(y,(m+1):2m)
    A_mul_B!(v, A.Q, u)
    return 0.0
end




function LinAlg.A_mul_B!{T,QType<:AbstractMatrix{T}}(Y::AbstractArray{T,1}, M::HSDEMatrix{T,QType}, B::AbstractArray{T,1})
    m2 = size(M)[1]
    mQ = size(M.Q)[1]
    @assert (m2,m2) == (size(Y,1), size(B,1))
    y1 = view(Y, 1:mQ)
    y2 = view(Y, (mQ+1):(2mQ))
    b1 = view(B, 1:mQ)
    b2 = view(B, (mQ+1):(2mQ))
    At_mul_B!(y1, M.Q, b2)
    A_mul_B!( y2, M.Q, b1)
    @inbounds y1 .= y1 .+ b1
    @inbounds y2 .= y2 .- b2
    return Y
end

LinAlg.At_mul_B!{T,QType<:AbstractMatrix{T}}(Y::AbstractArray{T,1}, Q::HSDEMatrix{T,QType}, B::AbstractArray{T,1}) =
    A_mul_B!(Y, Q, B)
