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
    debug::Int64
end

type NoStatus <: AbstractStatus
    status::Symbol
end

""" checkstatus(stat::Status, x)
Returns `false` if no check was made
If convergence check is done, returns `true` and sets stat.status to one of:
    :Continue, :Optimal, :Unbounded, :Infeasible
"""
function checkstatus(stat::Status, z; override = false)
    #TODO fix m, n
    if stat.i % stat.checki == 0 || override
        m, n, i, model, verbose, debug, ϵ = stat.m, stat.n, stat.i, stat.model, stat.verbose, stat.debug, stat.eps
        x, y, s, r, τ, κ = getvalues(model, m, n, z, i, verbose, debug)
        ϵpri = ϵdual = ϵgap = ϵinfeas = ϵunbound = ϵ
        p   = norm(model.A*x/τ + s/τ - model.b)/norm(1+norm(model.b))
        d   = norm(model.A'*y/τ + model.c - r/τ)/norm(1+norm(model.c))
        ctx = (model.c'*x)[1]
        bty = (model.b'*y)[1]
        g   = norm(ctx/τ + bty/τ)/(1+abs(ctx/τ)+abs(bty/τ))

        if debug > 0
            savedata(i, p, d, g, ctx, bty, κ, τ, x, y, s, model, debug)
        end
        # TODO test u^T v
        if verbose > 0
            printstatusiter(i, p, d, g, ctx, bty, κ/τ)
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

function checkstatus(stat::NoStatus, z)
    return false
end

function printstatusheader()
    println("-"^76)
    println(" Iter | pri res | dua res | rel gap | pri obj | dua obj | kap/tau")
    println("-"^76)
end

function printstatusiter(i, p, d, g, stx, bty, κτ)
    @printf("%6d|% .2e % .2e % .2e % .2e % .2e % .2e\n",i,p,d,g,stx,-bty,κτ)
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


@generated function savedata(i, p, d, g, ctx, bty, κ, τ, x, y, s, model, debug)
    ex = :(history = model.history)
    for v in  [:p, :d, :g, :ctx, :bty, :κ, :τ] #:t0
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
