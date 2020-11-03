using Convex
Random.seed!(2)

m = 40;  n = 50
A = randn(m, n); b = randn(m, 1)
x = Variable(n)
problem = minimize(sumsquares(A * x - b), [x >= 0])

ϵ = 1e-8
@static if VERSION >= v"1.5.0"
    # New random numbers in julia 1.5
    opt = 10.945929126466417
    #opt = 10.945829112949474
else
    opt = 12.38418747141913 # Before julia 1.5
end
aa = AndersonWrapper(DR(eps=ϵ, verbose=1))

solve!(problem, aa)

@test problem.status == :Optimal
@test problem.optval ≈ opt
@test abs(minimum(x.value)) < 10*ϵ

xsave = copy(x.value)

# Test indirect
problem = minimize(sumsquares(A * x - b), [x >= 0])
solve!(problem, GAPA(eps=1e-4, verbose=0))

@test problem.status == :Optimal
@test abs((problem.optval - opt)/opt) < 2e-3
@test maximum(abs.(x.value-xsave)) < 1e-3

#Test direct
problem = minimize(sumsquares(A * x - b), [x >= 0])
solve!(problem, GAPA(direct=true, eps=1e-4, verbose=0))

@test problem.status == :Optimal
@test abs((problem.optval - opt)/opt) < 2e-3
@test maximum(abs.(x.value-xsave)) < 1e-3

# Test β in GAPA
problem = minimize(sumsquares(A * x - b), [x >= 0])
solve!(problem, GAPA(0.5, 0.9, eps=1e-9, verbose=0))

@test problem.status == :Optimal
@test abs((problem.optval - opt)/opt) < 1e-8
@test maximum(abs.(x.value-xsave)) < 1e-7
