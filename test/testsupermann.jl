srand(2)

m = 480;  n = 500
A = sprandn(m, n,0.02); b = randn(m, 1)
x = Variable(n)
problem = minimize(sumsquares(A * x - b), [x >= 0])

ϵ = 1e-12
opt = 219.14056391716127
s = DR(direct=true, checki=100, max_iters=120000, eps=ϵ)
solve!(problem,s)

@test problem.status == :Optimal
@test problem.optval ≈ opt
@test abs(minimum(x.value)) < 10*ϵ

s = SuperMann(0., 0.99,0.99, 0.0005, 0.01,1.5, direct=true, checki=100, max_iters=10000, eps=ϵ)
solve!(problem,s)

@test problem.status == :Optimal
@test problem.optval ≈ opt rtol=10ϵ
@test abs(minimum(x.value)) < 10*ϵ
@test problem.model.enditr < 1000

s = SuperMann(direct=false, checki=100, max_iters=10000, eps=ϵ)
solve!(problem,s)

@test problem.status == :Optimal
@test problem.optval ≈ opt rtol=10ϵ
@test abs(minimum(x.value)) < 10*ϵ
@test problem.model.enditr < 1000
