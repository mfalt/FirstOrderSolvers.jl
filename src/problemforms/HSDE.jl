type HSDE{T,T2<:DualConeProduct}
    indAffine::IndAffine{T}
    indCones::T2
    n::Int64
end

function HSDE(model::FOSMathProgModel)
    m,n = size(model.A)
    Q = [spzeros(n,n) model.A'      model.c;
        -model.A      spzeros(m,m)  model.b;
        -model.c'     -model.b'     0      ]
    S1 = IndAffine([Q -speye(size(Q,1))], zeros(size(Q,1)))
    S2 = DualConeProduct(model.K1,model.K2)
    m,n = S2.m, S2.n
    status_generator = (mo, checki, eps, verbose, debug) ->
        Status(m, n, 0, mo, :Continue, checki, eps, verbose, false, debug)
    return HSDE{eltype(model.A), typeof(S2)}(S1, S2, 2*size(Q,1)), status_generator
end

function HSDE_getinitialvalue(model::FOSMathProgModel)
    m, n = size(model.A)
    l = m+n+1
    x = zeros(2*l)
    x[l]  = 1.0
    x[2l] = 1.0
    return x
end

function HSDE_populatesolution(model::FOSMathProgModel, x, status::Status)
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

HSDEStatus(model, checki, eps, verbose, debug) = Status(model.data.S2.m, model.data.S2.n, 0, model,
                                                        :Continue, checki, eps, verbose, false, debug)
