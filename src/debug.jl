function plothistory(problem, p = plot())
    hist = problem.model.history
    h = hist[:p]
    nonnan = !isnan.(hist[:p].values)
    val = hist[:p].values[nonnan]
    iter = hist[:p].iterations[nonnan]
    plot!(p, iter, val, yscale=:log10)
end
