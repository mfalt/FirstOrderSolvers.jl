mutable struct FeasibilityStatus <: AbstractStatus
    n::Int64
    i::Int64
    model::FeasibilityModel
    prev::Array{Float64,1} # Previous iterate to check convergence
    status::Symbol
    checki::Int64
    eps::Float64
    verbose::Int64
    checked::Bool
    direct::Bool
    init_time::UInt64     #As reported by time_ns()
    init_duration::UInt64 #In nano-seconds
    debug::Int64
end



""" checkstatus(stat::FeasibilityStatus, x)
Returns `false` if no check was made
If convergence check is done, returns `true` and sets stat.status to one of:
    :Continue, :Optimal, :Infeasible
"""
function checkstatus(stat::FeasibilityStatus, z; override = false)
    #TODO fix m, n
    if stat.i % stat.checki == 0 || override
        t = time_ns() - stat.init_time
        n, i, model, verbose, debug, ϵ = stat.n, stat.i, stat.model, stat.verbose, stat.debug, stat.eps
        prev = stat.prev
        # TODO Measure better
        err = norm(prev - z)
        if debug > 0
            savedata(i, err, z, t, model, debug)
        end
        if verbose > 0
            if !stat.direct
                cgiter = getcgiter(stat.model.data)
                push!(model.history, :cgiter, i, cgiter)
                printstatusiter(i, err, cgiter, t)
            else
                printstatusiter(i, err, t)
            end
        end

        status = :Continue
        if err <= ϵ
            if verbose > 0
                println("Found solution i=$i")
            end
            status = :Optimal
        elseif false # TODO Be able to check better
            status = :Infeasible
        end
        stat.status = status
        stat.checked = true
        stat.prev .= z
        return stat.checked
    else
        stat.checked = false
        stat.prev .= z
        return stat.checked
    end
end

function printstatusheader(stat::FeasibilityStatus)
    if stat.verbose > 0
        println("Time to initialize: $(stat.init_duration/1e9)s")
        width = 22 + (stat.direct ? 0 : 5)
        println("-"^width)
        print(" Iter | res")
        !stat.direct && print(" | cg ")
        println(" | time")
        println("-"^width)
    end
end

function printstatusiter(i, err, t)
    @printf("%6d|% 9.2e % .1es\n",i,err,t/1e9)
end

function printstatusiter(i, err, cgiter, t)
    @printf("%6d|% 9.2e % 4d % .1es\n",i,err,cgiter,t/1e9)
end

function savedata(i, err, z, t, model, debug)
    history = model.history
    push!(history, :err, i, err)
    push!(history, :t, i, t)
    if debug > 1
        push!(history, :z, i, z)
    end
    return
end
