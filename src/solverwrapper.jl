#TODO these are fallbacks
getinitialvalue(model::FOSMathProgModel, alg, data) =
    HSDE_getinitialvalue(model)

populate_solution(model::FOSMathProgModel, alg, data, x) =
    HSDE_populatesolution(model, x)


type Status
    m::Int64
    n::Int64
    i::Int64
    model::FOSMathProgModel
end

checkstatus(stat::Status, x) = stat.i % 100 == 0 && printstatus(stat.model, stat.m, stat.n, x, stat.i)

function solve!(model::FOSMathProgModel)
    opts = Dict(model.options)
    max_iters = :max_iters ∈ keys(opts) ? opts[:max_iters] : 10000
    x = getinitialvalue(model, model.alg, model.data)
    m,n = getmn(model.data) #TODO general
    status = Status(m,n,0,model)
    iterate(model.alg, model.data, status, x, max_iters)
    sol = populate_solution(model, model.alg, model.data, x)
    return sol
end


function iterate(alg::FOSAlgorithm, data::FOSSolverData, status, x, max_iters)
    t = tic()
    printstatusheader()
     #TODO general
    for i = 1:max_iters
        status.i = i
        step(alg, data, x, status)
    end
    println("Time for iterations: ")
    toc()
end


getmn(data::FOSSolverData) = data.S2.m, data.S2.n

function printstatus(model, m, n, z, i)
    nu = n+m+1
    x = view(z, 1:n)
    y = view(z, (n+1):(n+m))
    r = view(z, (nu+1):(nu+n))
    s = view(z, (nu+n+1):(nu+n+m))
    τ = z[nu]
    κ = z[2nu]
    printstatus(model, x, y, s, r, τ, κ, i)
end
