struct CGdata{T}
    r::Array{T,1}
    p::Array{T,1}
    z::Array{T,1}
    xinit::Array{T,1}#TODO remove this
    firstrun::Base.RefValue{Bool}
end

CGdata(r::AbstractArray{T},p::AbstractArray{T},z::AbstractArray{T}) where T = CGdata{T}(r,p,z, T[], Ref(true))

CGdata(size) = CGdata{Float64}(Array{Float64,1}(undef, size), Array{Float64,1}(undef, size), Array{Float64,1}(undef, size), Array{Float64,1}(undef, size), Ref(true))

function conjugategradient!(x,A,b; useasinitial=true)
    cgdata = CGdata(similar(b), similar(b), similar(b))
    if !useasinitial
        x .= 1.0
    end
    conjugategradient!(x, A, b, cgdata.r, cgdata.p, cgdata.z)
    return cgdata
end

"""
`conjugategradient!(x,A,b,r,p,Ap; tol = size(A,2)*eps(), max_iters = 10000)`

Solve `A*x==b` with the conjugate gradient method.
Uses `x` as warm start and stores solution in `x`.
`r`,`p`,`Ap` should be same size as `x` and will be changed.

Implemented as in "Matrix Compuatations" 2nd edition, Golub and Van Loan (1989)
"""
function conjugategradient!(x,A,b,r,p,Ap; tol = size(A,2)*eps(), max_iters = 10000)
    mul!(Ap, A, x)
    r .= b .- Ap
    p .= r
    rn = dot(r,r)
    iter = 1
    while true
        mul!(Ap, A, p)
        α = rn/dot(Ap,p)
        x .+= α.*p
        r .-= α.*Ap
        if norm(r) <= tol || iter >= max_iters
            break
        end
        rnold = rn
        rn = dot(r,r)
        β = rn/rnold
        #p .= r .+ β.*p
        p .*= β
        p .+= r
        iter += 1
    end
    iter == max_iters && @warn "CG reached max iterations, result may be inaccurate"
    return iter
end
