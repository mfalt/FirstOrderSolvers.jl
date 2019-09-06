#__precompile__()

module FirstOrderSolvers

using ProximalOperators
using SparseArrays
using Printf
import ValueHistories
using LinearAlgebra
import LinearAlgebra: mul!, Transpose, dot, norm

export Feasibility

include("cones.jl")
include("types.jl")
include("status.jl")
include("utilities/conjugategradients.jl")
include("utilities/affinepluslinear.jl")
include("problemforms/HSDE/HSDEStatus.jl")
include("problemforms/HSDE/HSDE.jl")
include("problemforms/HSDE/HSDEAffine.jl")

include("problemforms/Feasibility/Feasibility.jl")
include("problemforms/Feasibility/FeasibilityStatus.jl")

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
