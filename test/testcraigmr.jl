using Base.Test

include(joinpath(Pkg.dir("FirstOrderSolvers"), "src","utilities","CRAIGMR.jl"))

Base.norm(v::AbstractVector, M::AbstractMatrix) = sqrt(v'M*v)

function test_all_craig(M,N,A,m,n)
    τ = 1e-15; kmin = 10; kmax = 500
    K = [M A; A' -N]
    P = [M zeros(m,n);zeros(n,m) N]
    # M = isa(M, UniformScaling) ? M : factorize(M)
    # N = isa(N, UniformScaling) ? N : factorize(N)
    ### BEGIN TEST M != I N != I ########
    f = randn(m); g = zeros(n);
    xy = K\[f;g]
    x,y,i = craig_mr(A,M,N,f,kmin,τ,kmax)
    @test xy ≈ [x;y]
    @test norm(xy-[x;y],P) ≈ 0 atol = 1e-11 rtol=0


    #Test with g != 0, and some extra accuracy
    g = randn(n);
    xy = K\[f;g]
    x,y,i = craig_mr_full(A,M,N,f,g,kmin,τ,kmax)
    isave = i
    @test xy ≈ [x;y]
    println(norm(K*[x;y] -[f;g],P))
    @test norm(xy-[x;y],P) ≈ 0 atol = 1e-11 rtol=0
    ϵ = 1e-4
    #Test with g != 0 and with warmstart
    x0, y0 = x+ϵ*randn(m), y+ϵ*randn(n)
    x,y,i = craig_mr_warmstart(A,M,N,f,g,kmin,τ,kmax,x0,y0)
    @test xy ≈ [x;y]
    println(norm(K*[x;y] -[f;g],P))
    @test norm(xy-[x;y],P) ≈ 0 atol = 1e-11 rtol=0
    @test i < isave

    #Few iterations, good improvement
    x,y,i = craig_mr_warmstart(A,M,N,f,g,kmin,τ,20,x0,y0)
    println("$(norm(xy-[x;y],P)) <? $(norm(xy-[x0;y0],P))")
    println("$(norm(K*[x;y] -[f;g],P)) <? $(norm(K*[x0;y0] -[f;g],P))")
    @test norm(xy-[x;y],P) < norm(xy-[x0;y0],P)
    @test norm(K*[x;y] -[f;g],P) < norm(K*[x0;y0] -[f;g],P)
end

function test_all_craig_inplace(M,N,A,m,n)
    τ = 1e-15; kmin = 10; kmax = 500
    K = [M A; A' -N]
    P = [M zeros(m,n);zeros(n,m) N]
    # M = isa(M, UniformScaling) ? M : factorize(M)
    # N = isa(N, UniformScaling) ? N : factorize(N)
    ### BEGIN TEST M != I N != I ########
    f = randn(m); g = zeros(n);
    xy = K\[f;g]
    x, y = similar(f), similar(g)
    fsave = copy(f)
    i,craigdata = craig_mr!(x,y,A,M,N,f,kmin,τ,kmax)
    @test xy ≈ [x;y]
    @test norm(xy-[x;y],P) ≈ 0 atol = 1e-11 rtol=0
    @test f == fsave #Test that f was not overwritten

    #Test with g != 0, and some extra accuracy
    g = randn(n);
    gsave = copy(g)
    xy = K\[f;g]
    i,craigdata = craig_mr_full!(x,y,A,M,N,f,g,kmin,τ,kmax)
    isave = i
    @test xy ≈ [x;y]
    println(norm(K*[x;y] -[f;g],P))
    @test norm(xy-[x;y],P) ≈ 0 atol = 1e-11 rtol=0
    @test f == fsave #Test that f was not overwritten
    @test g == gsave #Test that g was not overwritten

    ϵ = 1e-4
    #Test with g != 0 and with warmstart
    y0 = copy(y)+ϵ*randn(n)
    x0 = copy(x)+ϵ*randn(m)
    x0save, y0save = copy(x0), copy(y0)
    i,craigdata = craig_mr_warmstart!(x,y,A,M,N,f,g,kmin,τ,kmax,x0,y0)
    @test xy ≈ [x;y]
    println(norm(K*[x;y] -[f;g],P))
    @test norm(xy-[x;y],P) ≈ 0 atol = 1e-11 rtol=0
    @test i < isave
    @test f == fsave    #Test that f was not overwritten
    @test g == gsave    #Test that g was not overwritten
    @test x0 == x0save  #Test that x0 was not overwritten
    @test y0 == y0save  #Test that y0 was not overwritten

    #Few iterations, good improvement
    i,craigdata = craig_mr_warmstart!(x,y,A,M,N,f,g,kmin,τ,20,x0,y0)
    println("$(norm(xy-[x;y],P)) <? $(norm(xy-[x0;y0],P))")
    println("$(norm(K*[x;y] -[f;g],P)) <? $(norm(K*[x0;y0] -[f;g],P))")
    @test norm(xy-[x;y],P) < norm(xy-[x0;y0],P)
    @test norm(K*[x;y] -[f;g],P) < norm(K*[x0;y0] -[f;g],P)
end

srand(10)
M = I; N = I;
m = 100; n = 200;
A = randn(m,n);
test_all_craig(M,N,A,m,n)
test_all_craig_inplace(M,N,A,m,n)

m = 200; n = 100;
A = randn(m,n);
test_all_craig(M,N,A,m,n)
test_all_craig_inplace(M,N,A,m,n)

M = I; N = I;
m = 200; n = 200;
A = randn(m,n);
test_all_craig(M,N,A,m,n)
test_all_craig_inplace(M,N,A,m,n)

m = 100; n = 200;
A = randn(m,n);
M1 = randn(m,m); M = sqrtm(M1*M1');
N1 = randn(n,n); N = sqrtm(N1*N1');
test_all_craig(M,N,A,m,n)
test_all_craig_inplace(M,N,A,m,n)

m = 200; n = 100;
A = randn(m,n);
M1 = randn(m,m); M = sqrtm(M1*M1');
N1 = randn(n,n); N = sqrtm(N1*N1');
test_all_craig(M,N,A,m,n)
test_all_craig_inplace(M,N,A,m,n)
