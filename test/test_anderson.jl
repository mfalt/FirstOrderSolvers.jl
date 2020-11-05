using ProximalOperators, Test, Random
Random.seed!(2)

xsol1 = randn(100)

A = randn(50,100)
b = A*xsol1

S1 = IndAffine(A,b)
S2 = IndBox(0.0, Inf)

problem = Feasibility(S1, S2, 100)

alg = DR(eps=1e-8, verbose=1, direct=true, checki=10)
sol1, model1 = FirstOrderSolvers.solve!(problem, alg, checki=10)
@test minimum(sol1.x) > -1e-12
@test maximum(abs, A*sol1.x - b) < 1e-12


alg = AndersonWrapper(DR(eps=1e-8, verbose=1, direct=true, checki=1))
sol2, model2 = FirstOrderSolvers.solve!(problem, alg, checki=1)
@test minimum(sol2.x) > -1e-12
@test maximum(abs, A*sol2.x - b) < 1e-12


using Convex, Random
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
aa = DR(eps=0.001*ϵ, direct=true)
solve!(problem, aa)

aa = AndersonWrapper(DR(eps=ϵ, direct=true))
solve!(problem, aa)

@test problem.status == :Optimal
@test problem.optval ≈ opt
@test abs(minimum(x.value)) < 10*ϵ

xsave = copy(x.value)

# Test indirect
problem = minimize(sumsquares(A * x - b), [x >= 0])
solve!(problem, GAPA(eps=1e-4, verbose=1))
solve!(problem, AndersonWrapper(GAPA(eps=1e-4, verbose=1)))

@test problem.status == :Optimal
@test abs((problem.optval - opt)/opt) < 2e-3
@test maximum(abs.(x.value-xsave)) < 1e-3

#Test direct
problem = minimize(sumsquares(A * x - b), [x >= 0])
solve!(problem, GAP(direct=true, eps=1e-4, verbose=1, max_iters=100000))
solve!(problem, AndersonWrapper(GAP(direct=true, eps=1e-4, verbose=1, max_iters=100000)))

@test problem.status == :Optimal
@test abs((problem.optval - opt)/opt) < 2e-3
@test maximum(abs.(x.value-xsave)) < 1e-3

# Test β in GAPA
problem = minimize(sumsquares(A * x - b), [x >= 0])
solve!(problem, GAPA(0.5, 0.9, eps=1e-9, verbose=0))

@test problem.status == :Optimal
@test abs((problem.optval - opt)/opt) < 1e-8
@test maximum(abs.(x.value-xsave)) < 1e-7
