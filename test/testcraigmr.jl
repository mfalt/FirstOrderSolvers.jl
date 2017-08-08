using Base.Test

include(joinpath(Pkg.dir("FirstOrderSolvers"), "src","utilities","CRAIGMR.jl"))

Base.norm(v::AbstractVector, M::AbstractMatrix) = sqrt(v'M*v)

function test_all_craig(M,N,A,m,n)
    τ = 1e-15; kmin = 10; kmax = 100
    K = [M A; A' -N]
    P = [M zeros(m,n);zeros(n,m) N]
    ### BEGIN TEST M != I N != I ########
    f = randn(m); g = zeros(n);
    xy = K\[f;g]
    x,y,i = craig_mr_simple(A,M,N,f,kmin,τ,5*kmax)
    @test xy ≈ [x;y]
    @test norm(xy-[x;y],P) ≈ 0 atol = 1e-12 rtol=0


    #Test with g != 0, and some extra accuracy
    g = randn(n);
    xy = K\[f;g]
    x,y,i = craig_mr_full(A,M,N,f,g,kmin,τ,5*kmax)
    @test xy ≈ [x;y]
    println(norm(K*[x;y] -[f;g],P))
    @test norm(xy-[x;y],P) ≈ 0 atol = 1e-12 rtol=0

    #Test with g != 0 and with warmstart
    x,y,i = craig_mr_warmstart(A,M,N,f,g,kmin,τ,2*kmax,x,y)
    @test xy ≈ [x;y]
    println(norm(K*[x;y] -[f;g],P))
    @test norm(xy-[x;y],P) ≈ 0 atol = 1e-12 rtol=0
    @test i < 50
end

srand(10)
M = I; N = I;
m = 100; n = 200;
A = randn(m,n);
test_all_craig(M,N,A,m,n)

m = 200; n = 100;
A = randn(m,n);
test_all_craig(M,N,A,m,n)

m = 100; n = 200;
A = randn(m,n);
M1 = randn(m,m); M = sqrtm(M1*M1');
N1 = randn(n,n); N = sqrtm(N1*N1');
test_all_craig(M,N,A,m,n)

m = 200; n = 100;
A = randn(m,n);
M1 = randn(m,m); M = sqrtm(M1*M1');
N1 = randn(n,n); N = sqrtm(N1*N1');
test_all_craig(M,N,A,m,n)
