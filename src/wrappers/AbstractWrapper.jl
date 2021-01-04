# Wrapper should contain field alg with original algorithm
abstract type AbstractWrapper <: FOSAlgorithm end

# WrapperData should contain field algdata with original algorithm data
abstract type AbstractWrapperData <: FOSSolverData end

# Default fallback for getting internal number of coordinate-gradient data
getcgiter(data::AbstractWrapperData) = getcgiter(data.algdata)

# Default fallback for getting solution
getsol(alg::AbstractWrapper, data::AbstractWrapperData, x) = getsol(alg.alg, data.algdata, x)

# Other functions needed for wrapper:

# function MyWrapper(alg::T; wrapperargs, wrapperkwargs..., kwargs...)
#     # Processing
#     MyWrapper(alg, wrapperkwargs..., merge(alg.options, kwargs))
# end

# init_algorithm!(ls::MyWrapper, model::AbstractFOSModel)

# Base.step(alg::MyWrapper, data::MyWrapperData, x, i, status::AbstractStatus, longstep=nothing)