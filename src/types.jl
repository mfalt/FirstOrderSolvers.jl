importall MathProgBase.SolverInterface
export FOSSolver

type Solution
    x::Array{Float64, 1}
    y::Array{Float64, 1}
    s::Array{Float64, 1}
    status::Symbol
end

abstract FOSSolverData


# Define Solver for interface
immutable FOSSolver <: AbstractMathProgSolver
    options
end
FOSSolver(;kwargs...) = FOSSolver(kwargs)


type FOSMathProgModel <: AbstractConicModel
    input_numconstr::Int64            # Only needed for interface?
    input_numvar::Int64               # Only needed for interface?
    K1::ConeProduct
    K2::ConeProduct
    A::SparseMatrixCSC{Float64,Int}   # The A matrix (equalities)
    b::Vector{Float64}                # RHS
    c::Vector{Float64}                # The objective coeffs (always min)
    alg::FOSAlgorithm                 #This descides which algorithm to call
    data::FOSSolverData               # Solver specific data is stored here
#    orig_sense::Symbol               # Original objective sense
    # Post-solve
    solve_stat::Symbol
    obj_val::Float64
    primal_sol::Vector{Float64}
    dual_sol::Vector{Float64}
    slack::Vector{Float64}
    options::Dict{Symbol, Any}
    enditr::Int64
    history::ValueHistories.MVHistory
end

function FOSMathProgModel(s::FOSAlgorithm; kwargs...)
    solverType = FOSMethods[FOSMethod]
    FOSMathProgModel{solverType}(0, 0, ConeProduct(), ConeProduct(), spzeros(0, 0),
                     Float64[], Float64[], solverType(), :NotSolved,
                     0.0, Float64[], Float64[], Float64[], Dict{Symbol,Any}(kwargs), -1, ValueHistories.MVHistory())
end

initialize!{T<:FOSMethodData}(model::FOSMathProgModel{T}) =
    error("No initialize! method defined for solver data type $T")
