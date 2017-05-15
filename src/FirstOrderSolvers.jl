__precompile__()

module FirstOrderSolver

using ProximalOperators
import ValueHistories

include("cones.jl")
include("types.jl")
include("FOSSolverInterface.jl")  # MathProgBase interface
include("solvers/GAP.jl")
include("solvers/GAPLS.jl")
include("solvers/GAPLSOld.jl")
include("solvers/GAPLSSplit.jl")
include("solvers/GAPSimple.jl")
include("solvers/SCS.jl")
include("lapack.jl")
include("solvermethods.jl")
include("line_search_methods.jl")

const FOSMethods = Dict{Symbol, Type}(
    :GAP => GAPData,
    :GAPLS => GAPLSDataSimple,
    :GAPLSOld => GAPLSDataOld,
    :GAPLSSplit => GAPLSSplitData,
    :GAPSimple => GAPDataSimple,
    :SCS => SCSData
)

end # module
