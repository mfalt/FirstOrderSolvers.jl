# Problematic problems from Convex.jl test files
using Convex

### Test 1
solver = DR(eps=1e-6, checki=100, max_iters=100000)
solver = GAPA(eps=1e-8, max_iters=10000)
solver = GAPA(eps=1e-6, max_iters=100000, checki=100, direct=true)
solver = LongstepWrapper(GAPA(eps=1e-6, max_iters=100000, checki=100, direct=true), longinterval=900, nsave=4)
solver = AP(eps=1e-8, max_iters=10000, checki=10, direct=true)

x = Variable(Positive())
y = Variable((3, 3))
p = minimize(x + y[1, 1], isposdef(y), x >= 1, y[2, 1] == 1)

x = Variable(Positive())

solve!(p, solver)
p.optval


solver = GAPA(eps=1e-8, max_iters=100000, checki=100, direct=true)
solver = DR(eps=1e-8, max_iters=100000, checki=10, direct=true)
solver = LongstepWrapper(GAPA(eps=1e-8, max_iters=10000, checki=10, direct=true), longinterval=50, nsave=2)
solver = LongstepWrapper(GAPA(eps=1e-8, max_iters=10000, checki=10, direct=true), longinterval=50, nsave=5)

using Random
x = Variable(200, 1)
Random.seed!(1)
A = randn(500,200)
b = randn(500)
p = minimize(norm2(A * x + b))
solve!(p, solver)
p.optval

v0 = p.optval
