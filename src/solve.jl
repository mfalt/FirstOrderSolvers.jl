
function solve!(model::FOSMathProgModel{GAPData}, K1, K2)
    t0 = time_ns()
    #Settings
    max_iters, ϵ, α, α₁, α₂, verbose, icheck, debug, adaptive = getOptions(model.options,
        [:max_iters, :eps, :alpha, :alpha1, :alpha2, :verbose, :icheck, :debug, :adaptive],
        (100_000,    1e-6, 0.8,    1.5,     1.5,     0,        20,      0,      100_000)
    )

    #Initialize variables
    m, n = size(model.A)
    nu = m+n+1
    nv = m+n+1
    nz = (m+n+1)*2

    z = model.data.z
    uμ = model.data.uμ
    z1 = model.data.z1
    z2 = model.data.z2
    # Initialize to avoid zero solutions
    ui = 1:nu
    vi = (nu+1):nz
    Q = model.data.Q
    F = model.data.F

    #z = [x;y;τ;r;s;κ]
    z1x = view(z1, 1:n)
    z1y = view(z1, (n+1):(n+m))
    z1r = view(z1, (nu+1):(nu+n))
    z1s = view(z1, (nu+n+1):(nu+n+m))
    z2x = view(z2, 1:n)
    z2y = view(z2, (n+1):(n+m))
    z2r = view(z2, (nu+1):(nu+n))
    z2s = view(z2, (nu+n+1):(nu+n+m))
    # K1 = model.K1
    # K2 = model.K2
    tmpcount = 0
    ϵpri = ϵdual = ϵgap = ϵinfeas = ϵunbound = ϵ
    primalGap = dualGap = linMove = coneMove = Inf
    status = :UnknownError
    if verbose > 0
        printstatusheader()
    end
    normNew = 1
    normOld = 1
    angle = .4
    for i = 1:max_iters
        #Update Affine Step ###########
        #uμ = F\z
        ldltsolve!(F, uμ, z)

        # z1[ui] = α₁ * uμ[ui] + (1-α₁)*z[ui]
        @inbounds @simd for j = 1:nu
            z1[j] = α₁ * uμ[j] + (1-α₁)*z[j]
        end

        #z1[vi] = Q*uμ[ui]*α₁ + (1-α₁)*z[vi]
        A_mul_B!(α₁, Q, view(uμ, ui), 0.0, view(z1, vi))
        @inbounds @simd for j in vi
            z1[j] += (1-α₁)*z[j]
        end

        # Update Cone step ############
        # Project C = K2×K1*×R+ ∈ R^n,R^m,R
        project!(    K2, z1x, z2x)
        projectDual!(K1, z1y, z2y)
        z2[nu] = max(z1[nu], 0)
        # Project C* = K2*×K1×R+ ∈ R^n,R^m,R
        projectDual!(K2, z1r, z2r)
        project!(    K1, z1s, z2s)
        z2[nz] = max(z1[nz], 0)

        if adaptive <= 2
            angNow = acos(dot(z1-z2,  z1-z)/(norm(z2-z1)*norm(z1-z)))
            angle = 0.*angle+1.*angNow
        end
        if adaptive <= 2# && angle < pi/100
            if i%100 == 0
                println(angle)
            end
            # α₁ = 2*(2-sqrt(2)*sqrt(1-cos(2*angle)))/(1+cos(2*angle))
            # α₂ = α₁
            #α₁ = 2
            #α₂ = 4*(1-sin(2*angle))/(1+cos(4*angle))
            α₁ = 2/(1+sin(angle))
            α₂ = α₁
            #α = 1.
            if adaptive == 1
                α = α₁/2#min(2*(sin(angle)+cos(angle))/(3*sin(angle)+cos(angle)), 1.0)
            end
            tmpcount = 0
        end
        # else
        #     if adaptive < max_iters && tmpcount > 50
        #         α₁ = 1.8
        #         α₂ = 1.8
        #         α = 1
        #     else
        #         tmpcount += 1
        #     end
        # end

        # Check before relaxing so cone constraints are satisfied if we quit
        if i%icheck == 0
            τ = z2[nu]
            κ = z2[nz]
            stat = checkStatus(model, z2x, z2y, z2s, z2r, τ, κ, ϵpri, ϵdual, ϵgap, ϵinfeas, ϵunbound, t0, i, verbose, debug)
            if stat != :Continue
                status = stat
                model.enditr = i
                break
            end
        end
        if i == max_iters
            status = :Indeterminate
            model.enditr = i
            break
        end

        # Relax α₂
        @inbounds @simd for j = 1:nz
            z2[j] = α₂*z2[j] + (1-α₂)*z1[j]
        end
        normOld = normNew
        normNew = norm(z2-z)
        #println("ang: $angle, α₂: $α₂, ratio: $(normNew/normOld)")
        #
        # if  normNew/normOld < 1. && i%50 == 0
        #     if normNew/normOld < 1.0*(α₁-1)
        #         α₁ = max(0.99*α₁, 1)
        #     else
        #         α₁ = 2/(1+sin(angle))
        #     end
        #     α₂ = α₁
        # end
        # if i%400 == 0
        #     α₁= 1.985
        #     α₂= 1.985
        # end

        if debug > 0
            push!(model.history, :normr, norm(z2-z))
        end

        # Relax α
        @inbounds @simd for j = 1:nz
            z[j]  = α*z2[j] + (1-α)*z[j]
        end


        # push!(model.history, :z, i, copy(z))
        # push!(model.history, :z1, i, copy(z1))
        # push!(model.history, :z2, i, copy(z2))
    end
    #TODO status
    τ = z2[nu]
    κ = z2[nz]
    solution = Solution(z2[1:n]/τ, z2[n+1:n+m]/τ, z2[nu+n+1:nu+n+m]/τ, status)
