import QPDAS
#import OSQP

#Type to store n linear equalities + 1 pair of linear equalities
mutable struct SavedPlanes{T}
    A::Array{T,2}
    b::Array{T,1}
    n::Int64
    neq::Int64
    nineq::Int64
end

function projectonnormals!(s::SavedPlanes{T},x,y) where T
    n = length(x)
    m = length(s.b)

    eqidx =               1:(s.neq*(s.n+1))
    ineqidx = (s.neq*(s.n+1)+1):(s.neq+s.nineq)*(s.n+1)

    A = s.A[  eqidx,:]
    b = s.b[eqidx]
    C = s.A[ineqidx,:]
    d = s.b[ineqidx]
    z = -x
    #
    qp = QPDAS.QuadraticProgram(BigFloat.(A), BigFloat.(b), BigFloat.(-C), BigFloat.(-d), BigFloat.(z), I, scaling=true, ϵ=1e-12)
    sol, val = QPDAS.solve!(qp)
    y .= sol


    # model = OSQP.Model()
    # M = [A;-C]
    # u = [b;-d]
    # l = [b;fill(-Inf, length(d))]
    # OSQP.setup!(model; P=SparseMatrixCSC{Float64}(I, n, n), l=l, A=sparse(M), u=u, verbose=false,
    #     eps_abs=0.001*eps(), eps_rel=0.001*eps(),
    #     eps_prim_inf=0.001*eps(), eps_dual_inf=0.001*eps(), max_iter=10000)
    #
    # OSQP.update!(model; q=z)
    # results = OSQP.solve!(model)
    # y .= results.x

    # println(size(s.A))
    # println("A:")
    # println(s.A[  eqidx,:])
    # println("b:")
    # println(s.b[  eqidx])
    # println("C:")
    # println(s.A[ineqidx,:])
    # println("d:")
    # println(s.b[ineqidx])
    # println("x:")
    # println(x)
    return false
end

function SavedPlanes(x::AbstractVector{T}, n::Int, neq, nineq) where T
    total = (n+1)*neq+(n+1)*nineq #Save all eq and last ineq
    SavedPlanes{T}(similar(x,total, length(x)), similar(x,total), n, neq, nineq)
end

#
# #Type to store n linear equalities + 1 pair of linear equalities
# mutable struct SavedPlanes{T}
#     A::Array{T,2}
#     b::Array{T,1}
#     n::Int64
#     neq::Int64
#     nineq::Int64
# end
#
# function projectonnormals!(s::SavedPlanes{T},x,y) where T
#     n = length(x)
#     m = length(s.b)
#     #println(s.A)
#     #println(s.b)
#     AAt = s.A*transpose(s.A)
#     AAtF = try
#         bkfact!(AAt)
#     catch
#         y .= x
#         #println("Non symmetric")
#         return true
#     end
#     tmp = copy(s.b)
#     # tmp = -Ax+b
#     LinearAlgebra.gemv!('N', -one(T), s.A, x, one(T), tmp)
#     # tmp = (AA')⁻¹*(-Ax+b)
#     try
#         ldiv!(AAtF,tmp)
#     catch
#         y .= x
#         println("Noninvertable")
#         return true
#     end
#     # y = A'(AA')⁻¹*(-Ax+b)
#     mul!(y, transpose(s.A), tmp)
#     # y = x + A'(AA')⁻¹*(-Ax+b) = (I-A'(AA')⁻¹*A)x + A'(AA')⁻¹b
#     y .= y .+ x
#     #println("in 2")
#     return false
# end
#
# function SavedPlanes(x::AbstractVector{T}, n::Int, neq, nineq) where T
#     total = (n+1)*neq+nineq #Save all eq and last ineq
#     SavedPlanes{T}(similar(x,total, length(x)), similar(x,total), n, neq, nineq)
# end

#Add planes to a specific location
function addplanesat(v1,b1,i,s::SavedPlanes)
    s.A[i,:] .= v1
    s.b[i]    = b1
end

#Add planes to a random location in 1:n with probability p=0.05
function addplanesrand(v1,b1,s::SavedPlanes, p = 0.05)
    if s.n > 0 && rand() < p
        i = rand(1:s.n)
        addplanesat(v1,b1,i,s)
    end
end
