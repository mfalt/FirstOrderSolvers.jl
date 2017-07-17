__precompile__()

module FirstOrderSolvers

using ProximalOperators
import ValueHistories

include("cones.jl")
include("types.jl")
include("status.jl")
include("utilities/conjugategradients.jl")
include("utilities/affinepluslinear.jl")
include("problemforms/HSDE.jl")
include("problemforms/HSDEAffine.jl")
include("FOSSolverInterface.jl")  # MathProgBase interface
include("solverwrapper.jl")

include("wrappers/longstep.jl")
include("wrappers/linesearch.jl")
include("solvers/defaults.jl")
include("solvers/solvers.jl")
#include("lapack.jl")
#include("solvermethods.jl")
#include("line_search_methods.jl")

end # module
