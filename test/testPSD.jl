using Convex, FirstOrderSolvers, SCS, Base.Test, ProximalOperators

ys = [-0.0064709 -0.22443;
     -0.22443    -1.02411]
y = Variable((2, 2))

#Solve with SCS
p = minimize(vecnorm(y-ys), isposdef(y))
solve!(p, SCSSolver(eps=1e-8, verbose=0))

ysol  = copy(y.value)

#Solve with projection
Y = Symmetric(ys)
X = similar(Y)
prox!(X, IndPSD(), Y)

#Testing SCS and projection equal
@test X.data ≈ ysol atol=1e-8

#Solve with DR
p = minimize(vecnorm(y-ys), isposdef(y))
solve!(p, DR(eps=1e-8, verbose=0))

@test y.value ≈ ysol atol=1e-8
