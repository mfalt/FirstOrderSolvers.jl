function printstatusheader()
    println("-"^76)
    println(" Iter | pri res | dua res | rel gap | pri obj | dua obj | kap/tau")
    println("-"^76)
end

function printstatusiter(i, p, d, g, stx, bty, κτ)
    @printf("%6d|% .2e % .2e % .2e % .2e % .2e % .2e\n",i,p,d,g,stx,-bty,κτ)
end


"""
Checks primal and dual feasibility assuming that x,y,s are satisfying the cone constraints
"""
function printstatus(model, x, y, s, r, τ, κ, i, verbose=1)
    p   = norm(model.A*x/τ + s/τ - model.b)/norm(1+norm(model.b))
    d   = norm(model.A'*y/τ + model.c - r/τ)/norm(1+norm(model.c))
    ctx = (model.c'*x)[1]
    bty = (model.b'*y)[1]
    g   = norm(ctx/τ + bty/τ)/(1+abs(ctx/τ)+abs(bty/τ))
    # TODO test u^T v
    if verbose > 0
        printstatusiter(i, p, d, g, ctx, bty, κ/τ)
    end
    return
end
