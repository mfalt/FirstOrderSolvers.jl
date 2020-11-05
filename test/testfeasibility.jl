using ProximalOperators
Random.seed!(2)

xsol1 = randn(100)

A = randn(50,100)
b = A*xsol1

S1 = IndAffine(A,b)
S2 = IndBox(0.0, Inf)

problem = Feasibility(S1, S2, 100)


sol, model = FirstOrderSolvers.solve!(problem, DR(eps=1e-8, verbose=1),checki=10)

@test sol.status == :Optimal
@test minimum(sol.x) > -1e-12
@test maximum(abs, A*sol.x - b) < 1e-12

sol, model = FirstOrderSolvers.solve!(problem, AP(eps=1e-8, verbose=0),checki=1)

@test sol.status == :Indeterminate

sol, model = FirstOrderSolvers.solve!(problem, GAP(eps=1e-8, verbose=0))

@test sol.status == :Indeterminate

sol, model = FirstOrderSolvers.solve!(problem, FISTA(eps=1e-8, verbose=0))

@test sol.status == :Indeterminate

@testset "Solver: $(typeof(alg))" for alg in
    [GAPP(eps=1e-8, verbose=0, proji=50),
    GAPA(eps=1e-8, verbose=0),
    LineSearchWrapper(GAP(eps=1e-8, verbose=0))]


    sol, model = FirstOrderSolvers.solve!(problem, alg)

    @test sol.status == :Optimal
    @test minimum(sol.x) > -1e-12
    @test maximum(abs, A*sol.x - b) < 1e-6
end

