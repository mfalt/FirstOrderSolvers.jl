"""
Representation of the matrix [I A'; A -I]
"""
immutable KKTMatrix{T, M<:AbstractMatrix{T}} <: AbstractMatrix{T}
    A::M
    am::Int64
    an::Int64
    # TODO test and remove if not better
    # x1::SubArray{T,1,Array{T,1},Tuple{UnitRange{Int64}},true}
    # x2::SubArray{T,1,Array{T,1},Tuple{UnitRange{Int64}},true}
    # y1::SubArray{T,1,Array{T,1},Tuple{UnitRange{Int64}},true}
    # y2::SubArray{T,1,Array{T,1},Tuple{UnitRange{Int64}},true}
end

KKTMatrix{T, M<:AbstractMatrix{T}}(A::M) = KKTMatrix{T,M}(A, size(A)...)#,
#            view([1.],1:1), view([1.],1:1), view([1.],1:1), view([1.],1:1))

# Alternative implementation with "constant" views

# @inbounds function Base.A_mul_B!{T}(y::AbstractVector{T}, M::KKTMatrix, x::AbstractVector{T})
#      if M.x1.parent !== x || M.y1.parent !== y
#         M.x1 = view(x,      1:M.an     )
#         M.x2 = view(x, M.an+1:M.an+M.am)
#         M.y1 = view(y,      1:M.an     )
#         M.y2 = view(y, M.an+1:M.an+M.am)
#         #println("Update pointers")
#     end
#
#     # [I A';
#     #  A -I]
#     At_mul_B!(M.y1, M.A, M.x2)
#     @blas! M.y1 += M.x1
#     A_mul_B!(M.y2,  M.A, M.x1)
#     @blas! M.y2 -= M.x2
# end

@inbounds function Base.A_mul_B!{T}(y::AbstractVector{T}, M::KKTMatrix, x::AbstractVector{T})
    x1 = view(x,      1:M.an     )
    x2 = view(x, M.an+1:M.an+M.am)
    y1 = view(y,      1:M.an     )
    y2 = view(y, M.an+1:M.an+M.am)

    # [I A';
    #  A -I]
    At_mul_B!(y1, M.A, x2)
    @blas! y1 += x1
    A_mul_B!(y2,  M.A, x1)
    @blas! y2 -= x2
end

#Assuming we don't introduce scaling ρ
Base.At_mul_B!{T}(y::AbstractVector{T}, M::KKTMatrix, x::AbstractVector{T}) = A_mul_B!(y, M, x)

"""
Representation of the function `f([x;z]) = q'x+i(Ax-βz==b)`
where `β` is `1` or `-1`.
"""
type AffinePlusLinear{T,M<:AbstractMatrix{T}} <: ProximalOperators.ProximableFunction
    M::KKTMatrix{T,M}
    A::M
    β::Int64
    b::Array{Float64,1}
    q::Array{Float64,1}
    rhs::Array{Float64,1}
    decreasing_accuracy::Bool
    i::Int64    #Count for number of calls, descides the tolerance
    cgiter::Int64 #How many iterarions CG did last time
    cgdata::CGdata{Float64}
end

function AffinePlusLinear{T,M<:AbstractMatrix{T}}(A::M, b, q, β; decreasing_accuracy=false)
    am, an = size(A)
    @assert β == 1 || β == -1
    #Initiate second half of rhs, maybe without view?
    rhs = Array{T}(am+an)
    rhs2 = view(rhs, (1+an):(am+an))
    rhs2 .= b
    return AffinePlusLinear{T,M}(KKTMatrix(A), A, β, b, q, rhs, decreasing_accuracy, 1, 0, CGdata(length(rhs)))
end

getcgiter(S::AffinePlusLinear) = S.cgiter

function prox!(y::AbstractArray, S::AffinePlusLinear, x::AbstractArray)
    #x2 = (2π-I)v or similar
    #[x;z] ∈ R(n×m)
    an = S.M.an
    am = S.M.am
    rhs1 = view(S.rhs, 1:an)
    x1   = view(x,     1:an)
    x2   = view(x,     (an+1):(an+am))
    y2   = view(y,     (an+1):(an+am))
    #Build rhs, rhs2 already contains b
    β = S.β
    At_mul_B!(rhs1, S.A, x2)
    rhs1 .= β.*rhs1 .+ x1 .- S.q
    #Now rhs = [x-q+β*A'z; b]

    #Do CG to solve M*[y1;y2] = rhs
    #Initial guess:
    cgdata = S.cgdata
    if cgdata.firstrun.x # has this been initialized? (pointer)
        cgdata.xinit .= x #Works as guess since M square
        cgdata.firstrun.x = false
    end
    #Since y is the initial guess, use previous solution
    y .= cgdata.xinit
    #Now run CG
    tol = if S.decreasing_accuracy
        max(0.2^sqrt(S.i), size(S.A,2)*eps())
    else
        size(S.A,2)*eps()
    end
    #println(tol)
    S.i += 1
    max_iters = 1000
    #TODO Verify that CG can solve this problem correctly
    iter = conjugategradient!(y, S.M, S.rhs, cgdata.r, cgdata.p, cgdata.z,
                                tol = tol, max_iters = max_iters)

    iter == max_iters && warn("CG reached max iterations, result may be inaccurate")
    S.cgiter = iter
    cgdata.xinit .= y #Save initial guess for next prox!
    #Now scale y2 with beta β
    y2 .*= β
    return 0.0
end
# A = [Ã -b;
#      0  1]
# [A -I][x;xn;z;zn] = 0
#
# [Ã -I][x;z] -b*xn = 0
# [Ã -I][x;z]  = b #OK
# xn-zn=0
#
# A' = [ Ã' 0;
#       -b' 1]
#
# hn=1?
# x̃n - xn -b'(z̃ - h) + 1(z̃n-hn) = 0
# [1 -b'][x̃n; z̃]=xn-b'h
# b'*x̃ =  b'*h #Maybe needed?
#
# #So:
# [I Ã'][x̃;z̃] = [-q+x+Ã'h]
# [b' 0][x̃;z̃] = [b'*h]
# #Complete system:
# [Ã -I ; = [b       ;
# [I  Ã'; =  -q+x+Ã'h;
# [b' 0 ] =  b'h     ]
