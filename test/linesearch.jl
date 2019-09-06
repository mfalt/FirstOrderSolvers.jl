using Convex
Random.seed!(2)

m = 40;  n = 50
A = randn(m, n); b = randn(m, 1)
x = Variable(n)
problem = minimize(sumsquares(A * x - b), [x >= 0])

ϵ = 1e-8
opt = 12.38418747141913
solver =  LineSearchWrapper(GAP(0.5, 1.0, 1.0, eps=ϵ, verbose=1, checki=1, max_iters=10000), lsinterval=100)
solve!(problem, solver)
