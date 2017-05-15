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
    S2 = DualConeProduct(K1,K2)
    HSDE{eltype(model.A), typeof(S2)}(S1, S2, 2*size(Q,1))
    return
end
