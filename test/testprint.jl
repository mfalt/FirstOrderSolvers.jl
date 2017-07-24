using Convex

function readmylines(rd)
    readline(rd) #Initialize time
    readline(rd) #Dashes
    l1 = readline(rd)
    readline(rd) #Dashes
    l2 = readline(rd)[1:7]
    l3 = readline(rd)[1:7]
    l4 = readline(rd)
    return l1, l2, l3, l4
end

#Expected outputs
o11 = " Iter | pri res | dua res | rel gap | pri obj | dua obj | kap/tau | cg  | time"
o12 = " Iter | pri res | dua res | rel gap | pri obj | dua obj | kap/tau | time"
o2  = "   100|"
o3  = "   200|"
o4  = "Found solution i=200"

srand(10)
n = 500
A = sprandn(n, 2n, 0.1)
x̄ = randn(2n)
b = A*x̄

x = Variable(2n)
p = minimize(norm(A*x-b), sum(x) == sum(x̄))

origstdout = STDOUT

#With cg output
s = GAPA(0.8, 0.9, direct=false, verbose=2, debug=0, eps=1e-8, checki=100)

rd,wr = redirect_stdout()
solve!(p,s)
l1, l2, l3, l4 = readmylines(rd)
redirect_stdout(origstdout)

@test l1 == o11
@test l2 == o2
@test l3 == o3
@test l4 == o4

@test evaluate(sum(x)-sum(x̄))       ≈ 0.0 atol=1e-8
@test evaluate(maximum(abs(A*x-b))) ≈ 0.0 atol=1e-8

#Without cg output
s = GAPA(0.8, 0.9, direct=true, verbose=2, debug=0, eps=1e-8, checki=100)

rd,wr = redirect_stdout()
solve!(p,s)
l1, l2, l3, l4 = readmylines(rd)
redirect_stdout(origstdout)

@test l1 == o12
@test l2 == o2
@test l3 == o3
@test l4 == o4

@test evaluate(sum(x)-sum(x̄))       ≈ 0.0 atol=1e-8
@test evaluate(maximum(abs(A*x-b))) ≈ 0.0 atol=1e-8
