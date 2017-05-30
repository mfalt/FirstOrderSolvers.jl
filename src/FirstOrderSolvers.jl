__precompile__()

module FirstOrderSolvers

using ProximalOperators
import ValueHistories

include("cones.jl")
include("types.jl")
include("problemforms/HSDE.jl")
include("FOSSolverInterface.jl")  # MathProgBase interface
include("solverwrapper.jl")
include("status.jl")
include("wrappers/longstep.jl")
include("wrappers/linesearch.jl")
include("solvers/gap.jl")
include("solvers/gapa.jl")
include("solvers/gapproj.jl")
#include("lapack.jl")
#include("solvermethods.jl")
#include("line_search_methods.jl")

end # module
