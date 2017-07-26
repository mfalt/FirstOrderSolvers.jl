__precompile__()

module FirstOrderSolvers

using ProximalOperators
import ValueHistories

#Types
include("cones.jl")
include("types.jl")
include("status.jl")

#Utils
include("utilities/conjugategradients.jl")
include("utilities/affinepluslinear.jl")
include("utilities/lbfgs.jl")

#Problem forms
include("problemforms/HSDE.jl")
include("problemforms/HSDEAffine.jl")


#Interface to Mathprogbase
include("FOSSolverInterface.jl")

#Main solver loop
include("solverwrapper.jl")

#Solver wrappers
include("wrappers/longstep.jl")
include("wrappers/linesearch.jl")

#Solvers
include("solvers/defaults.jl")
include("solvers/solvers.jl")

end # module
