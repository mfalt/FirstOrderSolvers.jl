__precompile__()

module FirstOrderSolvers
export GAP, GAPA, GAPP

using ProximalOperators
import ValueHistories

include("cones.jl")
include("types.jl")
include("status.jl")
include("problemforms/HSDE.jl")
include("FOSSolverInterface.jl")  # MathProgBase interface
include("solverwrapper.jl")

include("wrappers/longstep.jl")
include("wrappers/linesearch.jl")
include("solvers/solvers.jl")
#include("lapack.jl")
#include("solvermethods.jl")
#include("line_search_methods.jl")

end # module
