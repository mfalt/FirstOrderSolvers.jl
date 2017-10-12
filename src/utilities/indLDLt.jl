# struct IndGraphSparse2{T <: ProximalOperators.RealOrComplex, Ti} <: ProximableFunction
#   m::Int
#   n::Int
#   A::SparseMatrixCSC{T, Ti}
#   F::Base.SparseArrays.CHOLMOD.Factor{T} #LDL factorization
#
#   tmp::Array{T, 1}
#   tmpx::SubArray{T, 1, Array{T, 1}, Tuple{UnitRange{Int64}}, true}
# end
#
# function IndGraphSparse2{T,Ti}(A::SparseMatrixCSC{T,Ti})
#   m, n = size(A)
#   K = [speye(n) A'; A -speye(m)]
#
#   F = LinAlg.ldltfact(K)
#
#   tmp = Array{T,1}(m + n)
#   tmpx = view(tmp, 1:n)
#   tmpy = view(tmp, (n + 1):(n + m)) #second part is always zeros
#   fill!(tmpy, 0)
#
#   res = Array{T,1}(m + n)
#   return IndGraphSparse2{T,Ti}(m, n, A, F, tmp, tmpx)
# end
#
# function prox!{T}(
#     xy::AbstractVector{T},
#     f::IndGraphSparse2,
#     cd::AbstractVector{T}
#   )
#   #instead of res = [c + f.A' * d; zeros(f.m)]
#   At_mul_B!(f.tmpx, f.A, view(cd, (1+f.n):(f.n+f.m)))
#   f.tmpx .+= view(cd, 1:f.n)
#   # A_ldiv_B!(f.res, f.F, f.tmp) #is not working
#   xy .= f.F \ f.tmp #note here f.tmp which is m+n array
#   return 0.0
# end
using SparseLDLT

struct IndGraphSparse3{T <: ProximalOperators.RealOrComplex, Ti} <: ProximableFunction
  m::Int
  n::Int
  A::SparseMatrixCSC{T, Ti}
  F::SparseLDLT.LDLFactor#LDL factorization

  tmp::Array{T, 1}
  tmpx::SubArray{T, 1, Array{T, 1}, Tuple{UnitRange{Int64}}, true}
end

function IndGraphSparse3{T,Ti}(A::SparseMatrixCSC{T,Ti})
  m, n = size(A)
  K = [speye(n) A'; A -speye(m)]

  F = SparseLDLT.ldltfact(K)

  tmp = Array{T,1}(m + n)
  tmpx = view(tmp, 1:n)
  tmpy = view(tmp, (n + 1):(n + m)) #second part is always zeros
  fill!(tmpy, 0)

  res = Array{T,1}(m + n)
  return IndGraphSparse3{T,Ti}(m, n, A, F, tmp, tmpx)
end

function prox!{T}(
    xy::AbstractVector{T},
    f::IndGraphSparse3,
    cd::AbstractVector{T}
  )
  #instead of res = [c + f.A' * d; zeros(f.m)]
  At_mul_B!(f.tmpx, f.A, view(cd, (1+f.n):(f.n+f.m)))
  f.tmpx .+= view(cd, 1:f.n)
  # A_ldiv_B!(f.res, f.F, f.tmp) #is not working
  #xy .= f.F \ f.tmp #note here f.tmp which is m+n array
  A_ldiv_B!(xy, f.F, f.tmp)
  return 0.0
end
