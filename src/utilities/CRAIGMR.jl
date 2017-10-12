
""" Fifo queue that supports `sum(fun::Function, q::FIFO)` """
type FIFO{T}
    q::Array{T,1}
    n::Int64    #Total number of saved variables
    i::Int64    #Next index to be replaced in FIFO
end

""" Initialize fifo queue with NaN """
FIFO{I<:Integer}(n::I) = FIFO{Float64}(fill!(Array{Float64,1}(n), NaN), n, 1)

"""
    push!(q::FIFO{T}, val::T)

Push value `val` to queue `q`
"""
function Base.push!{T}(f::FIFO{T}, v::T)
    k = mod(f.i-1,f.n)+1
    f.q[k] = v
    f.i += 1
    return f
end

Base.sum(f::Union{Function, Type}, itr::FIFO) = sum(f, itr.q)
immutable CraigData
    m::Int64
    n::Int64
    d::Array{Float64,1}
    d̄::Array{Float64,1}
    rhs1::Array{Float64,1}
    rhs2::Array{Float64,1}
    u::Array{Float64,1}
    v::Array{Float64,1}
    tmp::Array{Float64,1}
    fifo::FIFO{Float64}
end

function CraigData(m,n,kmin)
    d = Array{Float64,1}(m)
    d̄ = Array{Float64,1}(m)
    rhs1 = Array{Float64,1}(m)
    rhs2 = Array{Float64,1}(n)
    u = Array{Float64,1}(m)
    v = Array{Float64,1}(n)
    tmp = Array{Float64,1}(n)   #This should only be used in craig_mr_full!
    fifo = FIFO(kmin)
    return CraigData(m,n,d,d̄,rhs1,rhs2,u,v,tmp,fifo)
end

"""
    x,y,k = craig_mr!(x,y,A,M,N,b,kmin,τ,kmax,craigdata::CraigData)

Solves problem `[M A;A'-N][x;y]=[b;0]` for `x` and `y` to accuracy `τ`.
Using at least `kmin` iterations and maxmum `kmax`.
"""
function craig_mr!(x,y,A,M,N,b,kmin,τ,kmax,craigdata::CraigData=CraigData(m,n,kmin))
    d, d̄, rhs1, rhs2, u, v = craigdata.d, craigdata.d̄, craigdata.rhs1,
                             craigdata.rhs2, craigdata.u, craigdata.v
    m,n = size(A)
    if m != craigdata.m || n != craigdata.n
        error("Invalid CraigData of size $((craigdata.m, craigdata.n)) supplied for system of size $((m,n))")
    end
    ζq = craigdata.fifo #FIFO Queue to store ζ values
    if kmin != ζq.n
        warn("Cragdata had invalid queue length for ζ, updating")
        ζq.q = Array{Float64,1}(kmin)
        ζq.n = kmin
        ζq.i = 1
    end
    fill!(ζq.q,NaN)

    fill!(d, 0.0)
    if all(iszero,b) #Fallback to avoid NaNs
        fill!(x, 0)
        fill!(y, 0)
        return 0, craigdata
    end
    #β1*M*u1 = b
    rhs1 .= b
    A_ldiv_B!(u,M,rhs1)     #u .= M\rhs1
    β = sqrt(u'rhs1)
    u ./= β
    #α1*N*v1 = A'u1
    At_mul_B!(rhs2, A, u)   #rhs2 .= A'u
    A_ldiv_B!(v,N,rhs2)     #v = N\rhs2
    α = sqrt(v'rhs2)
    v ./= α

    δ  = 1.;        αh = √(α^2+1);         c  = α/αh; s = 1/αh
    ζh = β;         αt = αh;               θ  = 0.
    d .= (1/αh).*u; d̄ .= 0.;               x .= 0.
    k  = 1;         Δ  = 0.;        converged = false
    while !converged && k < kmax
        #Start solve β{k+1}*M*u{k+1} = A*v{k} - α{k}*M*u{k}
        A_mul_B!(rhs1, M, u)    #rhs1 = M*u
        A_mul_B!(u, A, v)       #temporary store A*v in u
        rhs1 .= u .- α.*rhs1    #rhs1 .= A*v - α*M*u
        A_ldiv_B!(u,M,rhs1)     #u .= M\rhs1
        β = sqrt(u'rhs1)
        u ./= β                 # Now |u|_M = 1
        #End solve equation

        #Start solve α{k+1}*N*v{k+1} = A'u{k+1} - β{k+1}*N*v{k}
        A_mul_B!(rhs2, N, v)    # rhs2 = N*v
        At_mul_B!(v, A, u)      #temporary store A'u in v
        rhs2 .= v .- β.*rhs2    #rhs2 .= A'u - β*N*v
        A_ldiv_B!(v,N,rhs2)     #v = N\rhs2
        α = sqrt(v'rhs2)
        v ./= α                 # Now |v|_N = 1
        #End solve equation

        βh = c*β;          γ = s*β
        δ  = √(γ^2 +1);    c̄ = -1/δ; s̄ = γ/δ
        αh = √(α^2+δ^2);   c = α/αh; s = δ/αh
        ρ  = √(αt^2+βh^2); ĉ = αt/ρ; ŝ = βh/ρ
        d̄ .= (1/ρ).*(d .- θ.*d̄)               #Line 16 in book
        θ  = ŝ*αh;        αt = -ĉ*αh
        ζ  = ĉ*ζh;        ζh = ŝ*ζh; Δ = Δ + ζ^2
        d .= (1/αh).*(u .- βh.*d)
        #println("ζ^2: $(ζ^2), sum: $(sum(abs2,ζq)) , τ: $τ, τΔ: $(τ^2*Δ)")
        push!(ζq, ζ) #Push ζ to store last kmin values
        x .= x .+ ζ.*d̄
        if k ≥ kmin
            converged = sum(abs2,ζq) < τ^2*Δ
        end
        k = k + 1
    end
    # Calculate y .= N\(A'x)
    At_mul_B!(v,A,x) #Temporary use v
    A_ldiv_B!(y,N,v)
    return k,craigdata
end


"""
    x,y,k = craig_mr_full!(x,y,A,M,N,f,g,kmin,τ,kmax,craigdata=CraigData(m,n,kmin))

Solves problem `[M A;A'-N][x;y]=[f;g]` for `x` and `y`.

Overwrites `x` and `y`.
"""
function craig_mr_full!(x,y,A,M,N,f,g,kmin,τ,kmax,craigdata::CraigData=CraigData(m,n,kmin))
    y0, b = craigdata.tmp, craigdata.rhs1
    # x0 = 0
    #y0 = N\g
    A_ldiv_B!(y0, N, g)             # use tmp to store y0
    #b = f + A*y0
    A_mul_B!(b, A, y0) # use rhs1 to store b, this is safe
    b .+= f
    #println("f: $(norm(f)), b: $(norm(b))")
    #println("g: $(norm(g)), y0: $(norm(y0))")
    i,craigdata = craig_mr!(x,y,A,M,N,b,kmin,τ,kmax,craigdata)
    y .= y .- y0 # y = ŷ - y0
    return i,craigdata
end

"""
    x,y,k = craig_mr_warmstart!(x,y,A,M,N,f.g,kmin,τ,kmax,x0,y0)

Solves problem `[M A;A'-N][x;y]=[f;g]` for `x` and `y`.

Using initial guess x≈x0, y≈y0.

`x0` and `y0` can not point to same data as `x` and `y`.
"""
# function craig_mr_warmstart!(x,y,A,M,N,f,g,kmin,τ,kmax,x0,y0,craigdata::CraigData=CraigData(m,n,kmin))
#     f̂, ĝ = craigdata.rhs1, craigdata.rhs2
#     #f̂ = f - M*x0 - A*y0
#     A_mul_B!(x,A,y0)                   # use x to store A*y0
#     A_mul_B!(f̂, M, x0)                 # use rhs1 to store A*x0
#     f̂ .= f .- f̂ .- x
#     #ĝ = g - A'x0 + N*y0
#     A_mul_B!(y, N, y0)                 # use y to store N*y0
#     At_mul_B!(ĝ, A, x0)                # use rhs2 to store A'x0
#     ĝ .= g .- ĝ .+ y
#     τh = τ*norm(f)/norm(f̂) #Scale the tolerance by how good the guess is
#     i,craigdata = craig_mr_full!(x,y,A,M,N,f̂,ĝ,kmin,τh,kmax,craigdata) # xϵ, yϵ is same ref as x, y
#     x .= x0 .+ x # x = x0 + xϵ
#     y .= y0 .+ y # y = y0 + yϵ
#     return i,craigdata
# end
function craig_mr_warmstart!(x,y,A,M,N,f,g,kmin,τ,kmax,x0,y0,craigdata::CraigData=CraigData(m,n,kmin))
    f̂, ĝ = craigdata.u, craigdata.v
    #f̂ = f - M*x0 - A*y0
    A_mul_B!(x,A,y0)                   # use x to store A*y0
    A_mul_B!(f̂, M, x0)                 # use u to store A*x0
    f̂ .= f .- f̂ .- x
    #println("f: $(norm(f)), f̂: $(norm(f̂))")
    #ĝ = g - A'x0 + N*y0
    A_mul_B!(y, N, y0)                 # use y to store N*y0
    At_mul_B!(ĝ, A, x0)                # use v to store A'x0
    ĝ .= g .- ĝ .+ y
    #println("g: $(norm(g)), ĝ: $(norm(ĝ))")
    τh = τ*norm(f)/norm(f̂) #Scale the tolerance by how good the guess is
    i,craigdata = craig_mr_full!(x,y,A,M,N,f̂,ĝ,kmin,τh,kmax,craigdata) # xϵ, yϵ is same ref as x, y
    x .= x0 .+ x # x = x0 + xϵ
    y .= y0 .+ y # y = y0 + yϵ
    return i,craigdata
end

# Fallback for uniform scaling
function Base.A_ldiv_B!(y::AbstractVector, D::UniformScaling, x::AbstractVector)
    λinv = 1/D.λ
    y .= λinv .* x
end
function Base.A_mul_B!(y::AbstractVector, D::UniformScaling, x::AbstractVector)
    y .= D.λ .* x
end
function Base.At_mul_B!(y::AbstractVector, D::UniformScaling, x::AbstractVector)
    y .= D.λ .* x
end

# TODO TODO!!!
function Base.A_ldiv_B!(y::AbstractVector, D::AbstractMatrix, x::AbstractVector)
    y .= D\x
end
# function Base.A_mul_B!(y::AbstractVector, D::AbstractMatrix, x::AbstractVector)
#     y .= D.λ .* x
# end
# function Base.At_mul_B!(y::AbstractVector, D::AbstractMatrix, x::AbstractVector)
#     y .= D.λ .* x
# end


#
# """
#     x,y,k = craig_mr(A,M,N,b,kmin,τ,kmax)
#
# Solves problem `[M A;A'-N][x;y]=[b;0]` for `x` and `y` to accuracy `τ`.
# Using at least `kmin` iterations and maxmum `kmax`.
# """
# function craig_mr(A,M,N,b,kmin,τ,kmax)
#     m,n = size(A)
#     ζq = FIFO(kmin) #FIFO Queue to store ζ values
#
#     d = zeros(m)
#     #β1*M*u1 = b
#     rhs1 = copy(b)
#     u = M\rhs1
#     β = sqrt(u'rhs1)
#     u ./= β
#     #α1*N*v1 = A'u1
#     rhs2 = A'u
#     v = N\rhs2
#     α = sqrt(v'rhs2)
#     v ./= α
#
#     δ  = 1.;        αh = √(α^2+1);         c = α/αh; s = 1/αh
#     ζh = β;         αt = αh;               θ = 0.
#     d .= (1/αh).*u; d̄  = zeros(m);         x = zeros(m)
#     k  = 1;         Δ  = 0.;       converged = false
#     while !converged && k < kmax
#         #β{k+1}*M*u{k+1} = A*v{k} - α{k}*M*u{k}
#         rhs1 = A*v - α*M*u
#         u .= M\rhs1
#         β = sqrt(u'rhs1)
#         u ./= β
#         #α{k+1}*N*v{k+1} = A'u{k+1} - β{k+1}*N*v{k}
#         rhs2 = A'u - β*N*v
#         v = N\rhs2
#         α = sqrt(v'rhs2)
#         v ./= α
#
#         βh = c*β;          γ = s*β
#         δ  = √(γ^2 +1);    c̄ = -1/δ; s̄ = γ/δ
#         αh = √(α^2+δ^2);   c = α/αh; s = δ/αh
#         ρ  = √(αt^2+βh^2); ĉ = αt/ρ; ŝ = βh/ρ
#         d̄ .= (1/ρ).*(d .- θ.*d̄)               #Line 16 in book
#         θ  = ŝ*αh;        αt = -ĉ*αh
#         ζ  = ĉ*ζh;        ζh = ŝ*ζh; Δ = Δ + ζ^2
#         d .= (1/αh).*(u .- βh.*d)
#         push!(ζq, ζ) #Push ζ to store last kmin values
#         #println("ζ: $ζ, Δ: $Δ, sum: $(sum(abs2,ζq)) τ^2*Δ: $(τ^2*Δ)")
#         x .= x .+ ζ.*d̄
#         if k ≥ kmin
#             converged = sum(abs2,ζq) < τ^2*Δ
#         end
#         k = k + 1
#     end
#     y = N\(A'x)
#     return x,y,k
# end
#

"""
    x,y,k = craig_mr(A,M,N,b,kmin,τ,kmax)

Solves problem `[M A;A'-N][x;y]=[b;0]` for `x` and `y` to accuracy `τ`.
Using at least `kmin` iterations and maxmum `kmax`.

See also: `craig_mr!`
"""
function craig_mr(A,M,N,b,kmin,τ,kmax)
    x = similar(b)
    y = similar(b, size(A,2))
    i, craigdata = craig_mr!(x,y,A,M,N,b,kmin,τ,kmax)
    return x, y, i
end

# """
#     x,y,k = craig_mr_full(A,M,N,f,g,kmin,τ,kmax)
#
# Solves problem `[M A;A'-N][x;y]=[f;g]` for `x` and `y`.
# """
# function craig_mr_full(A,M,N,f,g,kmin,τ,kmax)
#     # x0 = 0
#     y0 = N\g
#     b = f + A*y0 # f - M*x0 - A*y0
#     x,ŷ,i = craig_mr(A,M,N,b,kmin,τ,kmax)
#     y = ŷ - y0
#     return x,y,i
# end


"""
    x,y,k = craig_mr_full(A,M,N,f,g,kmin,τ,kmax)

Solves problem `[M A;A'-N][x;y]=[f;g]` for `x` and `y`.

See also: `craig_mr_full!`
"""
function craig_mr_full(A,M,N,f,g,kmin,τ,kmax)
    x = similar(f)
    y = similar(g)
    i, craigdata = craig_mr_full!(x,y,A,M,N,f,g,kmin,τ,kmax)
    return x, y, i
end
"""
    x,y,k = craig_mr_warmstart(A,M,N,f,g,kmin,τ,kmax,x0,y0)

Solves problem `[M A;A'-N][x;y]=[f;g]` for `x` and `y`.

Using initial guess x≈x0, y≈y0.
"""
function craig_mr_warmstart(A,M,N,f,g,kmin,τ,kmax,x0,y0)
    f̂ = f - M*x0 - A*y0
    ĝ = g - A'x0 + N*y0
    τh = τ*norm(f)/norm(f̂) #Scale the tolerance by how good the guess is
    xϵ,yϵ,i = craig_mr_full(A,M,N,f̂,ĝ,kmin,τh,kmax)
    x = x0 + xϵ
    y = y0 + yϵ
    return x,y,i
end


# """
#     x,y,k = craig_mr_warmstart(A,M,N,f,g,kmin,τ,kmax,x0,y0)
#
# Solves problem `[M A;A'-N][x;y]=[f;g]` for `x` and `y`.
#
# Using initial guess x≈x0, y≈y0.
# """
# function craig_mr_warmstart(A,M,N,f,g,kmin,τ,kmax,x0,y0)
#     x = similar(f)
#     y = similar(g)
#     i, craigdata = craig_mr_warmstart!(x,y,A,M,N,f,g,kmin,τ,kmax,x0,y0)
#     return x, y, i
# end
