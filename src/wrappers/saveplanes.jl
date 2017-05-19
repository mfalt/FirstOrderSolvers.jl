
#Type to store n linear equalities + 1 pair of linear equalities
type SavedPlanes{T}
    A::Array{T,2}
    b::Array{T,1}
    n::Int64
    neq::Int64
    nineq::Int64
end

function projectonnormals!{T}(s::SavedPlanes{T},x,y)
    n = length(x)
    m = length(s.b)
    #println(s.A)
    #println(s.b)
    AAt = A_mul_Bt(s.A,s.A)
    AAtF = try
        bkfact!(AAt)
    catch
        y .= x
        #println("Non symmetric")
        return true
    end
    tmp = copy(s.b)
    # tmp = -Ax+b
    Base.LinAlg.gemv!('N', -one(T), s.A, x, one(T), tmp)
    # tmp = (AA')⁻¹*(-Ax+b)
    try
        A_ldiv_B!(AAtF,tmp)
    catch
        y .= x
        println("Noninvertable")
        return true
    end
    # y = A'(AA')⁻¹*(-Ax+b)
    At_mul_B!(y, s.A, tmp)
    # y = x + A'(AA')⁻¹*(-Ax+b) = (I-A'(AA')⁻¹*A)x + A'(AA')⁻¹b
    y .= y .+ x
    #println("in 2")
    return false
end

function SavedPlanes{T}(x::AbstractVector{T}, n::Int, neq, nineq)
    total = (n+1)*neq+nineq #Save all eq and last ineq
    SavedPlanes{T}(similar(x,total, length(x)), similar(x,total), n, neq, nineq)
end

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
