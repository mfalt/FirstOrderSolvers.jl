

function craig_mr_simple(A,M,N,b,kmin,τ,kmax)
    m,n = size(A)
    d = zeros(m)
    #β1*M*u1 = b
    rhs1 = copy(b)
    u = M\rhs1
    β = sqrt(u'u)#sqrt(u'rhs1)
    u ./= β
    #α1*N*v1 = A'u1
    rhs2 = A'u
    v = N\rhs2
    α = sqrt(v'v)#sqrt(v'rhs2)
    v ./= α

    δ  = 1.;        αh = √(α^2+1);         c = α/αh; s = 1/αh
    ζh = β;         αt = αh;               θ = 0.
    d .= (1/αh).*u; d̄  = zeros(m);         x = zeros(m)
    k  = 1;         Δ  = 0.;       converged = false
    while !converged && k < kmax
        #β{k+1}*M*u{k+1} = A*v{k} - α{k}*M*u{k}
        rhs1 = A*v - α*M*u
        u .= M\rhs1
        β = sqrt(u'u)#sqrt(u'rhs1)
        u ./= β
        #α{k+1}*N*v{k+1} = A'u{k+1} - β{k+1}*N*v{k}
        rhs2 = A'u - β*N*v
        v = N\rhs2
        α = sqrt(v'v)#sqrt(v'rhs2)
        v ./= α

        βh = c*β;          γ = s*β
        δ  = √(γ^2 +1);    c̄ = -1/δ; s̄ = γ/δ
        αh = √(α^2+δ^2);   c = α/αh; s = δ/αh
        ρ  = √(αt^2+βh^2); ĉ = αt/ρ; ŝ = βh/ρ
        d̄ .= (1/ρ).*(d .- θ.*d̄)               #Line 16 in book
        θ  = ŝ*αh;        αt = -ĉ*αh
        ζ  = ĉ*ζh;        ζh = ŝ*ζh; Δ = Δ + ζ^2
        d .= (1/αh).*(u .- βh.*d)

        x .= x .+ ζ.*d̄
        if k ≥ kmin
            converged = sum(abs2,ζ) < τ^2*Δ
        end
        k = k + 1
    end
    y = N\(A'x)
    return x,y
end
