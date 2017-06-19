using FirstOrderSolvers: conjugategradient!, CGdata

srand(2)
# a = randn(1000,1000)
# A = a'a
A = rand(1000,1000)
A = A'*A
b = randn(1000)
x0 = randn(1000)

x = copy(x0)
#anrm = sqrt(sum(abs2, A, 2))[:]
#A = A./anrm
#b = b./anrm
function test(A, b, x0, x, cgdata)
    x .= x0
    conjugategradient!(x, A, b, cgdata.r, cgdata.p, cgdata.z)
end

cgdata = CGdata(similar(b), similar(b), similar(b))
conjugategradient!(x, A, b, cgdata.r, cgdata.p, cgdata.z, max_iters = 100)
n0 = norm(A*x-b)    # 48
conjugategradient!(x, A, b, cgdata.r, cgdata.p, cgdata.z, max_iters = 5000)
n1 = norm(A*x-b)    # 3e-6

@test n1 < 1e-5

xcopy = x .+ 1e-5.*randn(size(x))
n2 = norm(A*xcopy-b) #0.9
conjugategradient!(xcopy, A, b, cgdata.r, cgdata.p, cgdata.z, max_iters = 100)
n3 = norm(A*xcopy-b) #3e-5

@test n3 < 10*n2
