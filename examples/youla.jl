using LinearAlgebra

N = 40
N_t = 300
N_f = 150
#s = tf("s");


""" `toeplitz2(vc, N = length(vc))`
Create toeplitxz matrix by stacking vc as columns, to width N, keeping height at length(vc) """
function toeplitz2(vc, N = size(vc,1))
    M = vc[:,:]
    Ncol = size(vc,2)
    Nrow = size(vc,1)
    for i = 2:N
        M = [M [zeros(i-1, Ncol); vc[1:Nrow-i+1, :]]]
    end
    return M
end

function toeplitz(vc::AbstractVector{T}, vr) where T
    N = length(vc)
    @assert N == length(vr)
    @assert vc[1] == vr[1]
    M = Array{T,2}(undef, N, N)
    # Subdiagonal
    for i = 1:N # Each column
        M[i:N,i] = vc[1:N-i+1]
    end
    for i = 1:N # Each row
        M[i, i:N] = vr[1:N-i+1]
    end
    return M
end

""" Compute A,b for the eqpression
    A*q+b corresponding to S*v=(I-PQ)v """
function Saff_t(TP,v,N)
    A = - T_P*toeplitz2(v,N)
    b = copy(v)
    return A, b
end

""" Compute A,b for the eqpression
    A*q+b corresponding to PS*v=P(I-PQ)v """
function PSaff_t(TP,v,N)
    AQ,bQ = CSaff_t(TP,v,N)
    A = -TP*TP*AQ
    b = TP*(v - TP*bQ)
    return A, b
end

""" Compute A,b for the eqpression
    A*q+b corresponding to CS*v=Qv """
function CSaff_t(TP,v,N)
    A = toeplitz2(v,N)
    b = zero(v)
    return A, b
end

""" Compute A,b for the eqpression
    A*q+b corresponding to CS*v=Qv """
function PCSaff_t(TP,v,N)
    A = TP*toeplitz2(v,N)
    b = zero(v)
    return A, b
end

function Saff_fr(P_fr, Qf)
    # Q_fr = G_q_fr*b
    # S_fr = 1 - P_fr.*Q_fr
    A = -P_fr.*Qf
    b = ones(eltype(Qf), size(P_fr,1))
    return A,b
end

""" T_P = SysP(g)
Given impulse response of system p, generate system (toeplitz) matrix TP """
function SysP(g)
    T_P = toeplitz(g, [g[1]; zeros(length(g)-1)]);
end

"""
    G_q_fr = Qfreq(Ω, N)
    Generate matrix for calculating frequency respose (Q_fr) of Q on grid Ω, such that
    Q_fr = G_q_fr*q
    where q is the impulse response or order N, and Ω∈(0,π)
"""
function Qfreq(Ω, N)
    exp.((-im.*Ω)*(0:N-1)')
end

using ControlSystems
#P = c2d(1/(20*s + 1) * exp(-5*s), 1);
z_d = tf("z", 1)
P_sysd = (1/20)/(z_d-0.95)/z_d^5
P_d(z) = (1/20)/(z-0.95)/z^5
#P(z) = (1/20)/(z-0.95)*z^(-5)

g = impulse(P_sysd, N_t-1)[1][:]

T_P = SysP(g)



Ω = 2*pi*(10 .^range(-3, log10(0.5), length=N_f))

#P_fr = freqresp(P_sysd, Ω)[:]
P_fr = P_d.(exp.(im.*Ω))
G_q_fr = Qfreq(Ω, N)

e1 = [1; zeros(N_t - 1)];
d1 = ones(N_t);

PSA, PSb = PSaff_t(T_P, d1, N)

SA_fr, Sb_fr = Saff_fr(P_fr, G_q_fr)

using Convex

################# Multiple affine and cones in Convex:

b = Variable(N)
y_d = Variable(N_t)
#Q_fr = ComplexVariable(N_f)
S_fr_re = Variable(N_f)
S_fr_im = Variable(N_f)
t1 = Variable()
t2 = Variable()
t3 = Variable(N_f)
t4 = Variable()
#c2 = (Q_fr == G_q_fr*b)

c3_re = (S_fr_re == real(SA_fr)*b + real(Sb_fr))        # Aff
c3_im = (S_fr_im == imag(SA_fr)*b + imag(Sb_fr))        # Aff

c4 = (t2 == t1 - 5^2)                                   # Aff
c5 = sumsquares(b) <= t1;                               # SOC    ∋ (b,t1)
c7 = (t2 <= 0)                                          # NonPos ∋ t2

c8 = (abs2(S_fr_re) + abs2(S_fr_im) <= t3)              # Many SOCs ∋ (S_fr,t3)
c9 = (t3 == 1.6^2)                                      # Affine

#cost = sumsquares(y_d)
c1 = (y_d == PSA*b+PSb)                                 # Affine
c10 = (sumsquares(y_d) <= t4)                           # SOC ∋ (y_d, t4)
c11 = (t4 == 1.1207+0.01)                               # Affine

prob = minimize(0, c1, c3_re, c3_im, c4, c5, c6, c7, c8, c9, c10, c11)

################# As one affine + cone in convex

Astag1 = toeplitz2([1 0 0; zeros(N_f-1, 3)])
Astag2 = toeplitz2([0 1 0; zeros(N_f-1, 3)])
Astag3 = toeplitz2([0 0 1; zeros(N_f-1, 3)])

A1 = [ zeros(N_f,1) -real(SA_fr)    zeros(N_f,1)    Astag2              zeros(N_f,1)    zeros(N_f, N_t) ]
A2 = [ zeros(N_f,1) -imag(SA_fr)    zeros(N_f,1)    Astag3              zeros(N_f,1)    zeros(N_f, N_t) ]
A3 = [   -1         zeros(1,N)      1               zeros(1,N_f*3)      zeros(1, 1)     zeros(1,   N_t) ]
A4 = [ zeros(N_f,1) zeros(N_f,N)    zeros(N_f,1)    Astag1              zeros(N_f,1)    zeros(N_f, N_t) ]
A5 = [ zeros(N_t,1) -PSA            zeros(N_t,1)    zeros(N_t,N_f*3)    zeros(N_t,1)    I               ]
A6 = [   0          zeros(1,N)      0               zeros(1,N_f*3)      1               zeros(1,   N_t) ]




x = [t1;b;t2;
    vcat(([t3[i];S_fr_re[i];S_fr_im[i];] for i in 1:N_f)...);
    t4; y_d]

# For comparison
xc = [sqrt(t1);b;t2;
    vcat(([sqrt(t3[i]);S_fr_re[i];S_fr_im[i];] for i in 1:N_f)...);
    sqrt(t4); y_d]

Aall = [A1;A2;A3;A4;A5;A6]
ball = [real(Sb_fr); imag(Sb_fr); -5^2; fill(1.6^2, N_f); PSb; 1.1207+0.01]

Aff = (Aall*x == ball)

Cone = [sumsquares(b) <= t1; t2 <= 0;
        vcat((abs2(S_fr_re[i]) + abs2(S_fr_im[i]) <= t3[i] for i in 1:N_f)...);
        sumsquares(y_d) <= t4]

prob = minimize(0, Aff, Cone...)
using SCS
@time solve!(prob, SCSSolver())
using Plots
plot(evaluate(b))



############ Using ProximalOperators
using ProximalOperators

ball2 = [real(Sb_fr); imag(Sb_fr); -5; fill(1.6, N_f); PSb; sqrt(1.1207+0.01)]

S1 = IndAffine(Aall, ball2)
fs = (IndSOC(), IndNonpositive(),
        (IndSOC() for i in 1:N_f)...,
        IndSOC())
idxs = (1:(N+1), (N+2):(N+2),
        ((N+3+i-1):(N+3+i+1) for i in 1:3:(N_f*3))...,
        (N+2+3*N_f+1):(N+2+3*N_f+N_t+1))
S2 = SlicedSeparableSum(fs, tuple.(idxs))

y0 = randn(size(Aall,2))
x0 = randn(size(Aall,2))

using FirstOrderSolvers
alg = DR()
problem = FirstOrderSolvers.Feasibility(S1, S2, size(Aall,2))
sol = FirstOrderSolvers.solve!(p, solver, checki=10, eps=1e-9)


using ProximalAlgorithms
solver = ProximalAlgorithms.DouglasRachford(gamma=1.0,maxit=10000,verbose=true)

x0, y0, it = solver(x0, f=S1, g=S2)
@profiler solver(x0, f=S1, g=S2)

tmp1 = similar(x0)
tmp2 = similar(x0)
tmp3 = similar(x0)
@time for i = 1:1000
    prox!(x0, S1, y0)
    tmp1 .= 2. *x0 .- y0
    prox!(tmp2, S2, tmp1)
    tmp3 .= 2. *tmp2 .- tmp1

    i%100 == 0 && (print("Err: "); println(norm(x0-tmp2)))
    i%100 == 0 && (print("Ang: "); println(acos(min(dot(y0-x0, tmp2-tmp1)/norm(x0-y0)/norm(tmp2-tmp1),1))))
    # Shadow in tmp2
    y0 .= (tmp3 .+ y0)./2
end

bsol = tmp2[2:N+1]
plot(bsol)

using Convex

b = Variable(N)

#T_Q = toeplitz2([b; zeros(N_t-N)])
#y_d = T_P*(Matrix{Float64}(I,N_t,N_t) - T_P*T_Q)*d1
#y_d = T_P*d1 - T_P*T_P*toeplitz2(d1,N)*b

y_d = PSA*b+PSb
Q_fr = G_q_fr*b
cost = sumsquares(y_d)
const1 = abs(SA_fr*b+Sb_fr) <= 1.6
const2 = sumsquares(b) <= 5^2
prob = minimize(cost, const1, const2)

using SCS
solve!(prob, SCSSolver())



using Plots
plot(evaluate(b))

plot(Ω./2 ./pi, abs.(1 .- P_fr.*evaluate(Q_fr)), yscale=:log10, xscale=:log10)

C_fr = vec(evaluate(Q_fr) ./ (1 .- P_fr .*evaluate(Q_fr)))

plot(Ω/2/pi, abs.(C_fr), yscale=:log10, xscale=:log10)
