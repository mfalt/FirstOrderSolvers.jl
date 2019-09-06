#Test against Convex.jl

using Convex

pkgdir(pkg::String) = abspath(joinpath(dirname(Base.find_package(pkg)), ".."))

single_solver_test_file = joinpath(pkgdir("Convex"),"test/runtests_single_solver.jl")
set_default_solver(DR(eps=1e-8, verbose=0, debug=0))

include(single_solver_test_file)

# TODO TODO
# #TODO Throws error, integrate with Base.Test
# FactCheck.exitstatus()
