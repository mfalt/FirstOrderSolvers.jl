
function normdiff(x::T,y::T) where T<:StridedVector
    nx, ny = length(x), length(y)
    s = 0.0
    @assert nx == ny
    for j = 1:nx
        s += abs2(x[j]-y[j])
    end
    return sqrt(s)
end