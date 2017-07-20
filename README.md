# FirstOrderSolvers

[![Build Status](https://travis-ci.org/mfalt/FirstOrderSolvers.jl.svg?branch=master)](https://travis-ci.org/mfalt/FirstOrderSolvers.jl)
[![Coverage Status](https://coveralls.io/repos/mfalt/FirstOrderSolvers.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/mfalt/FirstOrderSolvers.jl?branch=master)
[![codecov.io](http://codecov.io/github/mfalt/FirstOrderSolvers.jl/coverage.svg?branch=master)](http://codecov.io/github/mfalt/FirstOrderSolvers.jl?branch=master)

Package for large scale convex optimization solvers in julia. This package is intended to allow for easy **implementation**, **testing**, and **running** of solvers through the [Convex.jl](https://github.com/JuliaOpt/Convex.jl) interface.
The package is currently under **active development** and uses the [ProximalOperators.jl](https://github.com/kul-forbes/ProximalOperators.jl) package to do the low level projections.

## Installation
To run the solvers you need to have the following packages
```julia
Pkg.add("Convex")
Pkg.clone("https://github.com/mfalt/FirstOrderSolvers.jl.git")
```

## Usage
Define an optimization problem in the format supported by `Convex.jl`, and supply the desired solver to the `solve!` function. Exaple using DR for feasibility problems with the `GAP` solver 
```julia
using Convex, FirstOrderSolvers
m = 40;  n = 50
A = randn(m, n); b = randn(m, 1)
x = Variable(n)
problem = minimize(sumsquares(A * x - b), [x >= 0])

solve!(problem, GAP(0.5, 2.0, 2.0, max_iters=2000))
```

## Solvers
Currently, the available solvers are

| Solver | Description | Reference |
| --- | --- | --- |
| `GAP(α=0.8, α1=1.8, α2=1.8; kwargs...)` | Generalized Alternating Projections |    |
| `DR(α=0.5; kwargs...)` | Douglas-Rachford (`GAP(α, 2.0, 2.0)`)  | Douglas, Rachford (1956) |
| `AP(α=0.5; kwargs...)` | Alternating Projections (`GAP(α, 1.0, 1.0)`)  | Agmon (1954), Bregman (1967) |
| `GAPA(α=1.0; kwargs...)` | GAP Adaptive | [Fält, Giselsson (2017)](https://arxiv.org/abs/1703.10547) |
| `FISTA(α=1.0; kwargs...)` | FISTA |  Beck, Teboulle (2009) |
| `Dykstra(; kwargs...)` | Dykstra | Boyle, Dykstra (1986) |
| `GAPP(α=0.8, α1=1.8, α2=1.8; iproj=100; kwargs...)` | Projected GAP | [Fält, Giselsson (2016)](https://arxiv.org/abs/1609.05920)  |

## Keyword Arguments
All solvers accept for the following keyword arguments:

| Argument    | Default | Description (Values)                |
| ---         | ---     | ---                                 |
| `max_iters` | `10000` | Maximum number of iterations        |
| `verbose`   | `1`     |  Print verbosity level `0,1`        |
| `debug`     | `1`     | Level of debug data to save `0,1,2` |
| `eps`       | `1e-5`  | Accuracy of solution                |
| `checki`    | `100`   | Interval for checking convergence   |

## Debugging
If the keyword argument debug is set to `1` or `2` the following values will be stored in a [`ValueHistories.MVHistory`](https://github.com/JuliaPackageMirrors/ValueHistories.jl) in  `problem.model.history`, for each iteration the convergence check is run:

| Name        | Debug Level Required | Description |
| ---         | ---                  | ---                  |
| `:p`        | 1 | Relative Primal Residual |
| `:d`        | 1 | Relative Dual Residual |
| `:g`        | 1 | Relative Duality Gap |
| `:ctx`      | 1 | `cᵀx` |
| `:bty`      | 1 | `bᵀy` |
| `:κ`        | 1 | `κ` |
| `:τ`        | 1 | `τ` |
| `:x`        | 2 | `x` |
| `:y`        | 2 | `y` |
| `:s`        | 2 | `s` |

These values correspond to the values in the paper [Conic Optimization via Operator Splitting and Homogeneous Self-Dual Embedding (O'Donoghue et.al)](https://arxiv.org/abs/1312.3039).
