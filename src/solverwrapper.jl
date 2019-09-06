
function solve!(model::AbstractFOSModel)
    opts = Dict(model.options)
    #TODO better kwarg default handling
    max_iters = :max_iters ∈ keys(opts) ? opts[:max_iters] : 10000
    verbose = :verbose ∈ keys(opts) ? opts[:verbose] : 1
    debug = :debug ∈ keys(opts) ? opts[:debug] : 1
    eps = :eps ∈ keys(opts) ? opts[:eps] : 1e-5
    checki = :checki ∈ keys(opts) ? opts[:checki] : 100

    x = getinitialvalue(model, model.alg, model.data)
    #TODO general status
    status = model.status_generator(model, checki, eps, verbose, debug)
    guess = iterate(model.alg, model.data, status, x, max_iters)
    sol = populate_solution(model, model.alg, model.data, guess, status)
    return sol
end


function iterate(alg::FOSAlgorithm, data::FOSSolverData, status, x, max_iters)
    t1 = time()
    printstatusheader(status)
    for i = 1:max_iters
        status.i = i
        step(alg, data, x, i, status)
        if status.status != :Continue
            break
        end
    end
    #TODO Status should maybe be more consistent with guess
    guess = getsol(alg, data, x)
    if !status.checked #If max_iters reached without check
        checkstatus(status, guess, override = true) #Force check independent of iteration count
    end
    if status.verbose > 0
        println("Time for iterations: ")
        t2 = time()
        println("$(t2-t1) s")
    end
    return guess
end
