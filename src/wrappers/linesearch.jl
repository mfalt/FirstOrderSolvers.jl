export LineSearchWrapper

type LineSearchWrapper{T<:FOSAlgorithm} <: FOSAlgorithm
    lsinterval::Int64
    alg::T
    options
end

type LineSearchWrapperData{T<:FOSSolverData} <: FOSSolverData
    lsinterval::Int64
    tmp1::Array{Float64,1}
    tmp2::Array{Float64,1}
    tmp3::Array{Float64,1}
    mem::Array{Float64,1}
    res::Array{Float64,1}
    algdata::T
end

function LineSearchWrapper{T}(alg::T; lsinterval=100, kwargs...)
    if support_linesearch(alg) == Val{:False}
        error("Algorithm $T does not support line search")
    end
    LineSearchWrapper{T}(lsinterval, alg, [kwargs;alg.options])
end

function init_algorithm!(ls::LineSearchWrapper, model::FOSMathProgModel)
    alg, lsinterval = ls.alg, ls.lsinterval
    algdata, status =  init_algorithm!(alg, model)
    x = getinitialvalue(model, alg, algdata)
    data = LineSearchWrapperData(lsinterval, similar(x), similar(x), similar(x), similar(x), similar(x), algdata)
    return data, status
end

#getmn(data::LineSearchWrapperData) = getmn(data.algdata)

function Base.step(wrap::LineSearchWrapper, lsdata::LineSearchWrapperData, x, i, status::AbstractStatus, longstep=nothing)
    lsinterval, tmp1, tmp2, tmp3, mem, res =
        lsdata.lsinterval, lsdata.tmp1, lsdata.tmp2, lsdata.tmp3, lsdata.mem, lsdata.res

    do_linesearch = i % lsinterval == 0
    if do_linesearch  #Should we do ls?
        tmp1 .= x
        #step(wrap.alg, lsdata.algdata, x, i, status, nothing)

        #S1constant!(f0, wrap.alg, lsdata.algdata, mem, i)
        S1!(tmp2, wrap.alg, lsdata.algdata, x,    i, status, nothing)
        S2!(x,    wrap.alg, lsdata.algdata, tmp2, i, status, nothing)

        nostatus = NoStatus(:Continue)
        res .= x .- tmp1
        normres = norm(res)
        println("test, $normres")
        k = 1
        normresbest = Inf
        αbest = 1.0
        α = 0.1
        for k = 0:30
            α = α*1.8
            x .= tmp1 .+ α.*res
            #testres = get_normres!(wrap, lsdata, x, i, nostatus)
            S1!(tmp2, wrap.alg, lsdata.algdata, x,    i, nostatus, nothing)
            S2!(tmp3, wrap.alg, lsdata.algdata, tmp2, i, nostatus, nothing)
            testres = normdiff(x,tmp3)
            println("α: $α, $testres")
            if testres < normresbest
                normresbest = testres
                αbest = α
            end
        end
        println("α: $αbest")
        x .= tmp1 .+ αbest.*res
    else
        step(wrap.alg, lsdata.algdata, x, i, status, longstep)
    end
end

function normdiff{T<:StridedVector}(x::T,y::T)
    nx, ny = length(x), length(y)
    s = 0.0
    @assert nx == ny
    for j = 1:nx
        s += abs2(x[j]-y[j])
    end
    return sqrt(s)
end

function get_normres!(wrap::LineSearchWrapper, lsdata::LineSearchWrapperData, x, i, status::AbstractStatus)
    lsdata.tmp2 .= x
    step(wrap.alg, lsdata.algdata, x, i, status, nothing)
    #If time to do lsdata
    normres = normdiff(x,lsdata.tmp2)
end

function getsol(alg::LineSearchWrapper, data::LineSearchWrapperData, x)
    getsol(alg.alg, data.algdata, x)
end
