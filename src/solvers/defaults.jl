
support_longstep(::Any) = false
projections_per_step(::Any) = (0,0)

"""
Line search types:
Val{:False} : No line search available

Val{:True} : Line search available, but at the cost of extra evaluations
   Should provide an nonexpansive operator `y=T(x)` as
   `T!(y, alg<:FOSAlgorithm, data<:FOSSolverData, x, i, status, longstep)`
   so the algorithm becomes
   `x^(k+1) = (1-αk)x^k = αk*T(x^k)`
Val{:Fast} : Line search available with low extra cost
   Should provide two operators `y=S1(x)`, `z=S2(y)` as
   `S1!(y, alg<:FOSAlgorithm, data<:FOSSolverData, x, i, status, longstep)`
   `S2!(y, alg<:FOSAlgorithm, data<:FOSSolverData, x, i, status, longstep)`
   so the algorithm becomes
   `x^(k+1) = (1-αk)x^k = αk*S2(S1(x^k))`
   `S1` has to be an affine operator, i.e. S1(x+y)+S1(0)=S1(x)+S2(y).
"""
support_linesearch(::Any) = Val{:False}

#Defaults to that data.S1 counts cg iteration
function getcgiter(data::FOSSolverData)
    return getcgiter(data.S1)
end
