import ValueHistories: push!

type Status <: AbstractStatus
    m::Int64
    n::Int64
    i::Int64
    model::FOSMathProgModel
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


""" checkstatus(stat::Status, x)
Returns `false` if no check was made
If convergence check is done, returns `true` and sets stat.status to one of:
    :Continue, :Optimal, :Unbounded, :Infeasible
"""
function checkstatus(stat::Status, z; override = false)
    #TODO fix m, n
    if stat.i % stat.checki == 0 || override
        t = time_ns() - stat.init_time
        m, n, i, model, verbose, debug, ϵ = stat.m, stat.n, stat.i, stat.model, stat.verbose, stat.debug, stat.eps
        x, y, s, r, τ, κ = getvalues(model, m, n, z, i, verbose, debug)
        ϵpri = ϵdual = ϵgap = ϵinfeas = ϵunbound = ϵ
        p   = norm(model.A*x/τ + s/τ - model.b)/norm(1+norm(model.b))
        d   = norm(model.A'*y/τ + model.c - r/τ)/norm(1+norm(model.c))
        ctx = (model.c'*x)[1]
        bty = (model.b'*y)[1]
        g   = norm(ctx/τ + bty/τ)/(1+abs(ctx/τ)+abs(bty/τ))
        if debug > 0
            savedata(i, p, d, g, ctx, bty, κ, τ, x, y, s, t, model, debug)
        end
        # TODO test u^T v
        if verbose > 0
            if !stat.direct
                cgiter = getcgiter(stat.model.data)
                push!(model.history, :cgiter, i, cgiter)
                printstatusiter(i, p, d, g, ctx, bty, κ/τ, cgiter, t)
            else
                printstatusiter(i, p, d, g, ctx, bty, κ/τ, t)
            end
        end

        status = :Continue
        if p <= ϵpri*(1+norm(model.b)) && d <= ϵdual*(1+norm(model.c)) && g <= ϵgap*(1+norm(ctx/τ)+norm(bty/τ))
            if verbose > 0
                println("Found solution i=$i")
            end
            status =  :Optimal
        elseif norm(model.A*x + s) <= ϵunbound*(-ctx/norm(model.c))
            status = :Unbounded
        elseif norm(model.A'*y) <= ϵinfeas*(-bty/norm(model.b))
            status = :Infeasible
        end
        stat.status = status
        stat.checked = true
        return stat.checked
    else
        stat.checked = false
        return stat.checked
    end
end

function checkstatus(::NoStatus, z)
    return false
end

printstatusheader(::NoStatus) = nothing

function printstatusheader(stat::Status)
    if stat.verbose > 0
        println("Time to initialize: $(stat.init_duration/1e9)s")
        width = 76 + (stat.direct ? 0 : 5)
        println("-"^width)
        print(" Iter | pri res | dua res | rel gap | pri obj | dua obj | kap/tau")
        !stat.direct && print(" | cg ")
        println(" | time")
        println("-"^width)
    end
end

function printstatusiter(i, p, d, g, stx, bty, κτ, t)
    @printf("%6d|% 9.2e % 9.2e % 9.2e % 9.2e % 9.2e % 9.2e % .1es\n",i,p,d,g,stx,-bty,κτ,t/1e9)
end

function printstatusiter(i, p, d, g, stx, bty, κτ, cgiter, t)
    @printf("%6d|% 9.2e % 9.2e % 9.2e % 9.2e % 9.2e % 9.2e % 4d % .1es\n",i,p,d,g,stx,-bty,κτ,cgiter,t/1e9)
end

function getvalues(model, m, n, z, i, verbose, debug)
    nu = n+m+1
    x = view(z, 1:n)
    y = view(z, (n+1):(n+m))
    r = view(z, (nu+1):(nu+n))
    s = view(z, (nu+n+1):(nu+n+m))
    τ = z[nu]
    κ = z[2nu]
    x, y, s, r, τ, κ
end


# """
# Checks primal and dual feasibility assuming that x,y,s are satisfying the cone constraints
# """
# function printstatus(model, x, y, s, r, τ, κ, i, verbose, debug)
#     p   = norm(model.A*x/τ + s/τ - model.b)/norm(1+norm(model.b))
#     d   = norm(model.A'*y/τ + model.c - r/τ)/norm(1+norm(model.c))
#     ctx = (model.c'*x)[1]
#     bty = (model.b'*y)[1]
#     g   = norm(ctx/τ + bty/τ)/(1+abs(ctx/τ)+abs(bty/τ))
#     if debug > 0
#         savedata(i, p, d, g, ctx, bty, κ, τ, x, y, s, model, debug)
#     end
#     # TODO test u^T v
#     if verbose > 0
#         printstatusiter(i, p, d, g, ctx, bty, κ/τ)
#     end
#     return
# end


@generated function savedata(i, p, d, g, ctx, bty, κ, τ, x, y, s, t, model, debug)
    ex = :(history = model.history)
    for v in  [:p, :d, :g, :ctx, :bty, :κ, :τ, :t]
        ex = :($ex ; push!(history, $(parse(":$v")), i, $v))
    end
    ex2 = :()
    for v in [:x,:y,:s]
        ex2 = :($ex2 ; push!(history, $(parse(":$v")), i, copy($v)))
    end
    ex = :($ex; if debug > 1; $ex2; end)
    ex = :($ex; return)
    return ex
end
