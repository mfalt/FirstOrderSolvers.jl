# FirstOrderSolvers

[![Build Status](https://travis-ci.org/mfalt/FirstOrderSolvers.jl.svg?branch=master)](https://travis-ci.org/mfalt/FirstOrderSolvers.jl)
[![Coverage Status](https://coveralls.io/repos/mfalt/FirstOrderSolvers.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/mfalt/FirstOrderSolvers.jl?branch=master)
[![codecov.io](http://codecov.io/github/mfalt/FirstOrderSolvers.jl/coverage.svg?branch=master)](http://codecov.io/github/mfalt/FirstOrderSolvers.jl?branch=master)

Package for large scale convex optimization solvers. This package is intended to allow for easy implementation, testing, and running or solvers through the [Convex.jl](https://github.com/JuliaOpt/Convex.jl) interface.
The package is currently under active development.

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

solve!(problem, GAP(0.5,2.0,2.0, max_iters=2000))
```
Currently, the only available solvers are:
Generalized Alternating Projections:
```julia
GAP(α=0.8, α1=1.8, α2=1.8; kwargs...)
```
GAP Adaptive
```julia
GAPA(α=1.0; kwargs...)
```
Projected GAP
```julia
GAPP(α=0.8, α1=1.8, α2=1.8; iproj=100, kwargs...)
```

