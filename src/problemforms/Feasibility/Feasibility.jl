## Problem
struct Feasibility{T1<:ProximableFunction,T2<:ProximableFunction}
    S1::T1
    S2::T2
    n::Int64
end

mutable struct FeasibilitySolution <: AbstractSolution
    x::Array{Float64, 1}
    status::Symbol
end


# Model with all data
mutable struct FeasibilityModel{T1<:ProximableFunction,T2<:ProximableFunction} <: FOSModel
    S1::T1
    S2::T2
    n::Int64
    alg::FOSAlgorithm                 #This descides which algorithm to call
    data::FOSSolverData               # Solver specific data is stored here
    solve_stat::Symbol
    obj_val::Float64
    options::Dict{Symbol, Any}
    enditr::Int64
    status_generator::Function
    init_duration::UInt64           # In nano-seconds
    history::ValueHistories.MVHistory
end

function FeasibilityModel(problem::Feasibility, alg::FOSAlgorithm, kwargs...)
    t1 = time_ns()

    model = FeasibilityModel(problem.S1, problem.S2, problem.n, alg, FOSSolverDataPlaceholder(), :NotSolved,
                     0.0, Dict{Symbol,Any}(kwargs), -1,
                     x -> (@error "No status generator defined"),
                     UInt64(1), ValueHistories.MVHistory())

    data, status_generator = init_algorithm!(model.alg, model)
    model.status_generator = status_generator
    model.data = data
    t2 = time_ns()
    model.init_duration = t2-t1

    return model
end

function solve!(problem::Feasibility, alg::FOSAlgorithm; kwargs...)
    model = FeasibilityModel(problem, alg, kwargs...)
    solution = solve!(model)
end

getinitialvalue(model::FeasibilityModel) = zeros(model.n)
getinitialvalue(model::FeasibilityModel, alg, data) = zeros(model.n)

function populate_solution(model::FeasibilityModel, alg, data, x, status)
    endstatus = status.status
    if endstatus == :Continue
        endstatus = :Indeterminate
    end
    solution = FeasibilitySolution(x, endstatus)
end



"""
    `S1, S2, n, status_generator = get_sets_and_status(alg::FOSAlgorithm, model::FOSMathProgModel)`
 Get the sets S1, S2 and their size n
"""
function get_sets_and_status(alg::FOSAlgorithm, model::FeasibilityModel)
    direct = true # alg.direct # No support in this model for setting flag here
    #!alg.direct && @warn "Not possible to set direct/indirect solve in algorithm, use appropriate set from ProximalOperators"
    status_generator = (mo, checki, eps, verbose, debug) ->
        FeasibilityStatus(mo.n, 0, mo, fill(NaN, mo.n), :Continue, checki, eps, verbose, false, direct, time_ns(), mo.init_duration, debug)
    return model.S1, model.S2, model.n, status_generator
end
