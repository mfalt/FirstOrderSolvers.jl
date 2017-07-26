using FirstOrderSolvers: LBFGSstate, set_x_∇f!, H∇f!

srand(1)
m = 1000
Q = randn(m,m)
Q = (I+Q'Q)/4000
c = randn(m)
d = randn()

f(x) = (x-c)'Q*(x-c) + d
g(x) = Q*x - (Q*c)

const Qc = Q*c
function g!(y,x)
    A_mul_B!(y, Q, x)
    y .-= Qc
end
x0 = randn(m)

#x1save = Array{Float64,2}(m, 1000)
x1 = copy(x0)
for i = 1:1000
    x1 -= g(x1)
    #x1save[:,i] .= x1
end

#plot([f(x1save[:,i]) for i = 1:1000])

@test f(x1) ≈ f(c) atol = 0.02

#With H0 = I
n = 30
lbfgs = LBFGSstate(m,30, usegamma=false)

#x2save = Array{Float64,2}(m, 300)
x2 = copy(x0)
grad = copy(x2)
#grad = g(x2)
#initlbfgs!(x2, grad, lbfgs)
#x2 -= grad
for i = 1:300
    #println(norm(x2))
    g!(grad, x2)
    set_x_∇f!(x2, grad, lbfgs)  #Update state
    H∇f!(grad, lbfgs)           #Get newton direction
    x2 -= grad                  # x2 = x2 - H*∇f(x2)
    #x2save[:,i] .= x2
end

#plot!([f(x2save[:,i]) for i = 1:300])
@test f(x2) ≈ f(c) atol = 0.02
@test norm(x2-c)/norm(x2) < 0.01

#With H0 = γI
n = 30
lbfgs = LBFGSstate(m,30, usegamma=false)

#x2save = Array{Float64,2}(m, 300)
x2 = copy(x0)
grad = copy(x2)
#grad = g(x2)
#initlbfgs!(x2, grad, lbfgs)
#x2 -= grad
for i = 1:300
    #println(norm(x2))
    g!(grad, x2)
    set_x_∇f!(x2, grad, lbfgs)  #Update state
    H∇f!(grad, lbfgs)           #Get newton direction
    x2 -= grad                  # x2 = x2 - H*∇f(x2)
    #x2save[:,i] .= x2
end

#plot!([f(x2save[:,i]) for i = 1:300])
@test f(x2) ≈ f(c) atol = 0.02
@test norm(x2-c)/norm(x2) < 0.01
