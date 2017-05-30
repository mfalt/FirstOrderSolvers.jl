#TODO these are fallbacks
getinitialvalue(model::FOSMathProgModel, alg, data) =
    HSDE_getinitialvalue(model)

populate_solution(model::FOSMathProgModel, alg, data, x) =
    HSDE_populatesolution(model, x)


type Status <: AbstractStatus
    m::Int64
    n::Int64
    i::Int64
    model::FOSMathProgModel
    status::Symbol
    checki::Int64
    eps::Float64
    verbose::Int64
    debug::Int64
end

type NoStatus <: AbstractStatus
status::Symbol
end

function solve!(model::FOSMathProgModel)
    opts = Dict(model.options)
    max_iters = :max_iters ∈ keys(opts) ? opts[:max_iters] : 10000
    verbose = :verbose ∈ keys(opts) ? opts[:verbose] : 1
    debug = :debug ∈ keys(opts) ? opts[:debug] : 1
    eps = :eps ∈ keys(opts) ? opts[:eps] : 1e-5
    checki = :checki ∈ keys(opts) ? opts[:checki] : 1
    x = getinitialvalue(model, model.alg, model.data)
    m,n = getmn(model.data) #TODO general
    status = Status(m, n, 0, model, :Continue, checki, eps, verbose, debug)
    guess = iterate(model.alg, model.data, status, x, max_iters)
    sol = populate_solution(model, model.alg, model.data, guess)
    return sol
end


function iterate(alg::FOSAlgorithm, data::FOSSolverData, status, x, max_iters)
    t = tic()
    printstatusheader()
     #TODO general
    for i = 1:max_iters
        status.i = i
        step(alg, data, x, i, status)
        if status.status != :Continue
            break
        end
    end
    guess = getsol(alg, data, x)
    println("Time for iterations: ")
    toc()
    return guess
end


getmn(data::FOSSolverData) = data.S2.m, data.S2.n
