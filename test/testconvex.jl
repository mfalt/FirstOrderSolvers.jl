#Test against Convex.jl

using Convex

single_solver_test_file = joinpath(Pkg.dir("Convex"),"test/runtests_single_solver.jl")
set_default_solver(DR(eps=1e-8, verbose=0, debug=0))

include(single_solver_test_file)

# TODO TODO
# #TODO Throws error, integrate with Base.Test
# FactCheck.exitstatus()
