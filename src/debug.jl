function plothistory(hist, p = plot())
    h = problem.model.history[:p]
    nonnan = !isnan.(problem.model.history[:p].values)
    val = problem.model.history[:p].values[nonnan]
    iter = problem.model.history[:p].iterations[nonnan]
    plot!(p, iter, val, yscale=:log10)
end
