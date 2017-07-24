export LongstepWrapper

include("saveplanes.jl")

type LongstepWrapper{T<:FOSAlgorithm} <: FOSAlgorithm
    longinterval::Int64
    nsave::Int64
    alg::T
    options
end

type LongstepWrapperData{T<:FOSSolverData} <: FOSSolverData
    longinterval::Int64
    nsave::Int64
    saved::SavedPlanes{Float64}
    saveineq::Bool
    savepos::Int64
    eqi::Int64
    uneqi::Int64
    tmp::Array{Float64,1}
    algdata::T
end

function LongstepWrapper{T}(alg::T; longinterval=100, nsave=10, kwargs...)
    LongstepWrapper{T}(longinterval, nsave, alg, [kwargs;alg.options])
end

function init_algorithm!(long::LongstepWrapper, model::FOSMathProgModel)
    alg, longinterval, nsave = long.alg, long.longinterval, long.nsave
    !support_longstep(alg) && error("Algorithm alg does not support longstep")
    data, status_generator =  init_algorithm!(alg, model)
    neq, nineq = projections_per_step(alg)
    #TODO x
    x = HSDE_getinitialvalue(model)
    saved = SavedPlanes(x, nsave, neq, nineq)
    data = LongstepWrapperData(longinterval, nsave, saved, false, 0, 1, 1, similar(x), data)
    return data, status_generator
end

getmn(data::LongstepWrapperData) = getmn(data.algdata)

function Base.step(wrap::LongstepWrapper, longstep::LongstepWrapperData, x, i, status::AbstractStatus)
    nsave, longinterval = longstep.nsave, longstep.longinterval
    #println("At iteration $i")
    savepos = (i-1)%longinterval - longinterval + nsave + 2 # Index of collection to save
    #println("savepos: $savepos")
    if 0 < savepos  #Should we start saving collections?
        longstep.savepos = savepos
        longstep.eqi, longstep.uneqi = 0, 0
    end
    step(wrap.alg, longstep.algdata, x, i, status, longstep)
    #If time to do longstep
    if longstep.savepos == nsave+1
        #println("Doing projection!")
        fail = projectonnormals!(longstep.saved, x, longstep.tmp)
        longstep.savepos = -1
        x .= longstep.tmp
    end
end

function getsol(alg::LongstepWrapper, data::LongstepWrapperData, x)
    getsol(alg.alg, data.algdata, x)
end

addprojeq(::Void, ::Any, ::Any) = nothing
addprojineq(::Void, ::Any, ::Any) = nothing

function addprojeq{T<:LongstepWrapperData}(long::T, y, x)
    if long.savepos > 0
        #println("adding EQ, nsave: $(long.nsave) i:$(long.i) eqi:$(long.eqi) uneqi:$(long.uneqi)")
        s = long.saved
        i = (long.savepos-1)*(s.neq + s.nineq*long.saveineq) + long.eqi + 1
        s.A[i,:] .= x .- y
        b = 0.0
        for j = 1:length(x)
            b += (x[j]-y[j])*y[j]
        end
        s.b[i] = b
        long.eqi += 1
    end
    return
end

"""Saves inequalities as equality for now"""
function addprojineq{T<:LongstepWrapperData}(long::T, y, x)
    if long.savepos > 0 && (long.saveineq || long.savepos == long.nsave+1)
        #println("adding inEQ, nsave: $(long.nsave) i:$(long.i) eqi:$(long.eqi) uneqi:$(long.uneqi)")
        s = long.saved
        i = (long.savepos-1)*(s.neq + s.nineq*long.saveineq) + long.eqi + long.uneqi + 1
        s.A[i,:] .= x .- y
        b = 0.0
        for j = 1:length(x)
            b += (x[j]-y[j])*y[j]
        end
        s.b[i] = b
        long.uneqi += 1
    end
end
