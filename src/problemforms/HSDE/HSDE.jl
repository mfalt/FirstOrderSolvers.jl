mutable struct HSDE{T1<:ProximableFunction,T2<:DualConeProduct}
    indAffine::T1
    indCones::T2
    n::Int64
end

function HSDE(model::FOSMathProgModel; direct=false)
    m,n = size(model.A)

    S1 = if direct
            # Uses direct solver from ProximalOperators
            Q = [spzeros(n,n) model.A'      model.c;
                -model.A      spzeros(m,m)  model.b;
                -model.c'     -model.b'     0      ]
            IndAffine([sparse(Q) -I], zeros(size(Q,1)))
        else
            Q = HSDEMatrixQ(model.A, model.b, model.c)
            # Q = [spzeros(n,n) model.A'      model.c;
            #     -model.A      spzeros(m,m)  model.b;
            #     -model.c'     -model.b'     0      ]
            # Using CG on KKT system
            AffinePlusLinear(Q, zeros(size(Q,1)), zeros(size(Q,1)), 1, decreasing_accuracy=true)
        end
    S2 = DualConeProduct(model.K1,model.K2)
    m,n = S2.m, S2.n
    status_generator = (mo, checki, eps, verbose, debug) ->
        HSDEStatus(m, n, 0, mo, :Continue, checki, eps, verbose, false, direct, time_ns(), model.init_duration, debug)
    return HSDE(S1, S2, 2*(size(model.A,1)+size(model.A,2)+1)), status_generator
end

# TODO Unclear why this is nessesary/differenet
getinitialvalue(model::FOSMathProgModel) = HSDE_getinitialvalue(model)

getinitialvalue(model::FOSMathProgModel, alg, data) =
    HSDE_getinitialvalue(model)

populate_solution(model::FOSMathProgModel, alg, data, x, status) =
    HSDE_populatesolution(model, x, status)

function HSDE_getinitialvalue(model::FOSMathProgModel)
    m, n = size(model.A)
    l = m+n+1
    x = zeros(2*l)
    x[l]  = 1.0
    x[2l] = 1.0
    return x
end

function HSDE_populatesolution(model::FOSMathProgModel, x, status::HSDEStatus)
    m, n = size(model.A)
    l = m+n+1
    @assert length(x) == 2l
    τ = x[l]
    κ = x[2l]
    #TODO Eveluate status
    endstatus = status.status
    if endstatus == :Continue
        endstatus = :Indeterminate
    end
    solution = Solution(x[1:n]/τ, x[n+1:n+m]/τ, x[l+n+1:l+n+m]/τ, endstatus)
end
