importall MathProgBase.SolverInterface

ConicModel(s::FOSAlgorithm) = FOSMathProgModel(s; s.options...)
LinearQuadraticModel(s::FOSAlgorithm) = ConicToLPQPBridge(ConicModel(s))

function optimize!(m::FOSMathProgModel)
    #TODO fix history
    m.history = ValueHistories.MVHistory()
    solution = solve!(m) #TODO Code here!

    m.solve_stat = solution.status
    m.primal_sol = solution.x

    m.dual_sol = solution.y

    m.slack = solution.s

    m.obj_val = dot(m.c, m.primal_sol)# * (m.orig_sense == :Max ? -1 : +1)
end

status(m::FOSMathProgModel) = m.solve_stat
getobjval(m::FOSMathProgModel) = m.obj_val
getsolution(m::FOSMathProgModel) = copy(m.primal_sol)

function loadproblem!(model::FOSMathProgModel, c, A, b, constr_cones, var_cones)
    loadproblem!(model, c, sparse(A), b, constr_cones, var_cones)
end

function loadproblem!(model::FOSMathProgModel, c, A::SparseMatrixCSC, b, constr_cones, var_cones)
    t1 = time_ns()
    model.input_numconstr = size(A,1)
    model.input_numvar = size(A,2)

    # Verify only good cones
    for cone_vars in constr_cones
        cone_vars[1] in badcones && error("Cone type $(cone_vars[1]) not supported")
    end
    for cone_vars in var_cones
        cone_vars[1] in badcones && error("Cone type $(cone_vars[1]) not supported")
    end

    conesK1 = tuple([conemap[t[1]] for t in constr_cones]...)
    indexK1 = tuple([t[2] for t in constr_cones]...)
    model.K1 = ConeProduct(indexK1, conesK1)

    conesK2 = tuple([conemap[t[1]] for t in var_cones]...)
    indexK2 = tuple([t[2] for t in var_cones]...)
    model.K2 = ConeProduct(indexK2, conesK2)

    model.A = A
    # TODO figure out :Min/:Max
    model.b = b
    model.c = c

    # Calls a specific method based on the type T in model::FOSMathProgModel{T}
    data, status_generator = init_algorithm!(model.alg, model)
    model.status_generator = status_generator
    model.data = data
    t2 = time_ns()
    model.init_duration = t2-t1
    return model
end

numvar(model::FOSMathProgModel) = model.input_numvar
numconstr(model::FOSMathProgModel) = model.input_numconstr

supportedcones(s::FOSAlgorithm) = [:Free, :Zero, :NonNeg, :NonPos, :SOC, :SDP, :ExpPrimal, :ExpDual]
