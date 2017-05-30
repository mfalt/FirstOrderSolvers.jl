type LineSearchWrapper{T<:FOSAlgorithm} <: FOSAlgorithm
    lsinterval::Int64
    alg::T
    options
end

type LineSearchWrapperData{T<:FOSSolverData} <: FOSSolverData
    lsinterval::Int64
    tmp1::Array{Float64,1}
    tmp2::Array{Float64,1}
    res::Array{Float64,1}
    algdata::T
end

function LineSearchWrapper{T}(alg::T; lsinterval=100, kwargs...)
    LineSearchWrapper{T}(lsinterval, alg, kwargs)
end

function init_algorithm!(ls::LineSearchWrapper, model::FOSMathProgModel)
    alg, lsinterval = ls.alg, ls.lsinterval
    #!support_lsstep(alg) && error("Algorithm alg does not support ls")
    data =  init_algorithm!(alg, model)
    #TODO x
    x = HSDE_getinitialvalue(model)
    LineSearchWrapperData(lsinterval, 0, similar(x), similar(x), similar(x), data)
end

getmn(data::LineSearchWrapperData) = getmn(data.algdata)

function Base.step(wrap::LineSearchWrapper, lsdata::LineSearchWrapperData, x, i, status::AbstractStatus)
    lsinterval, tmp1, tmp2, res = lsdata.lsinterval, lsdata.tmp1, lsdata.tmp2, lsdata.res

    if i % lsinterval == 0  #Should we do ls?
        tmp1 .= x
    end
    step(wrap.alg, lsdata.algdata, x, i, status)

    nostatus = NoStatus(:Continue)
    if i % lsinterval == 0
        res .= x .- tmp1
        normres = norm(res)
        k = 1
        normresbest = normres
        αbest = 1.0
        for k = 1:20
            α = 2.0^k
            x .= tmp1 .+ α.*res
            testres = get_normres!(wrap, lsdata, x, nostatus)
            if testres < normresbest
                normresbest = testres
                αbest = α
            end
        end
        println("α: $αbest")
        x .= tmp1 .+ αbest.*res
    end
end

function get_normres!(wrap::LineSearchWrapper, lsdata::LineSearchWrapperData, x, status)
    lsdata.tmp2 .= x
    step(wrap.alg, lsdata.algdata, x, status)
    #If time to do lsdata
    normres = norm(x - lsdata.tmp2)
end

function getsol(alg::LineSearchWrapper, data::LineSearchWrapperData, x)
    getsol(alg.alg, data.algdata, x)
end
