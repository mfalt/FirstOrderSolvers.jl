using FirstOrderSolvers: KKTMatrix, AffinePlusLinear
using Convex

srand(10)

#Test KKTMatrix
A = randn(10,20)
M1 = [I A';A -I]
M2 = KKTMatrix(A)
x = randn(30)
y1 = randn(30)
y2 = randn(30)
A_mul_B!(y1, M1, x)
A_mul_B!(y2, M2, x)
@test y1 ≈ y2

At_mul_B!(y1, M1, x)
At_mul_B!(y2, M2, x)
@test y1 ≈ y2

#Test AffinePlusLinear

x0 = randn(20)
z0 = randn(10)
q = randn(20)
b = randn(10)

# β = 1
β = 1
S2 = AffinePlusLinear(A, b, q, β)
y2 = Array{Float64,1}(30)
FirstOrderSolvers.prox!(y2, S2, [x0;z0])
x2 = view(y2,1:20)
z2 = view(y2,21:30)

# #This test is equivalent to y3
# x = Variable(20)
# z = Variable(10)
# p = minimize(1/2*norm(x-x0)^2+1/2*norm(z-z0)^2+dot(q,x), [A*x-z==b])
# solve!(p, FirstOrderSolvers.DR(eps=1e-9))
# x1 = x.value
# z1 = z.value
# y1 = [x1; z1]
# @test y1 ≈ y2

y3 = [I A'; A -I]\[x0-q+A'z0; b]
@test y3 ≈ y2


β = -1
S2 = AffinePlusLinear(A, b, q, β)
y2 = Array{Float64,1}(30)
FirstOrderSolvers.prox!(y2, S2, [x0;z0])
x2 = view(y2,1:20)
z2 = view(y2,21:30)

# #This test is equivalent to y3
# x = Variable(20)
# z = Variable(10)
# p = minimize(1/2*norm(x-x0)^2+1/2*norm(z-z0)^2+dot(q,x), [A*x-β*z==b])
# solve!(p, FirstOrderSolvers.DR(eps=1e-9))
# x1 = x.value
# z1 = z.value
# y1 = [x1; z1]
# @test y1 ≈ y2

y3 = [I -A'; A I]\[x0-q-A'z0; b]
@test y3 ≈ y2
