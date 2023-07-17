x(t) = sinh(sinh(t))
dx_dt(t) = cosh(sinh(t)) * cosh(t)

function energy_output_integrals(sol, i, Xs, θ, par, H, H2)
    N = size(Xs, 1)

    Z = exp(sol[end, 1] - θ * H(view(Xs, :, i), par)) * √abs(extract_det_jac(sol.u[1], N))

    output = [Z, Z * H(sol[N+1:2N, 1], par)]

    regularize(output), false
end

function energy_integrals(θ, par, f!, H, lb, ub; 
    reltol=1e-2, abstol=1e-3, dif_eq_alg=nothing, int_alg=CubatureJLh(), 
    maxiters = 10^6, callback = caustic_callback)
    function formated_integrand(Xs, _)
        #@show size(Xs)
        integrand(Xs, θ, par, f!, H, energy_output_integrals;
        dif_eq_alg, reltol, abstol, save_start=false, save_everystep=false, verbose=false, callback)
    end

    prob = IntegralProblem(formated_integrand, lb, ub; nout=2, batch=Threads.nthreads())

    U = solve(prob, int_alg; reltol, abstol, maxiters)
    U[2] / U[1]
end

function energy_integrals2(θ, par, f!, H, lb, ub; 
    reltol=1e-2, abstol=1e-3, dif_eq_alg=nothing, int_alg=CubatureJLh(), 
    maxiters = 10^6, callback = caustic_callback)

    function formated_integrand(t, _)
        #@show size(t)
        jacs = prod(dx_dt,t,dims=1) |> vec
        map!(x,t,t)
        integ = integrand(t, θ, par, f!, H, energy_output_integrals;
        dif_eq_alg, reltol, abstol, save_start=false, save_everystep=false, verbose=false, callback)

        @tullio integ[m,n] *= jacs[n]

        integ
    end

    prob = IntegralProblem(formated_integrand, lb, ub; nout=2, batch=Threads.nthreads())

    U = solve(prob, int_alg; reltol, abstol, maxiters)
    U[2] / U[1]
end


function heat_output_integrals(sol, i, Xs, θ, par, H, H2)
    N = size(Xs, 1)

    Z = exp(sol[end, 1] - θ * H(view(Xs, :, i), par)) * √abs(extract_det_jac(sol.u[1], N))

    output = [Z, Z * H(sol[N+1:2N, 1], par), Z * H2(sol[N+1:2N, 1], par)]

    regularize(output), false
end

function heat_integrals(θ, par, f!, H, lb, ub;  
    reltol=1e-2, abstol=1e-3, dif_eq_alg=nothing, int_alg=CubatureJLh(), 
    maxiters = 10^6, callback = caustic_callback)
    formated_integrand(Xs, _) = integrand(Xs, θ, par, f!, H, heat_output_integrals;
        dif_eq_alg, reltol, abstol, save_start=false, save_everystep=false, verbose=false, callback)

    #pxs = [[0.1,0.1,x,y] for y ∈ LinRange(lb[4],ub[4],128), x ∈ LinRange(lb[3],ub[3],128)] |> vec |> stack
    #formated_integrand(pxs,nothing)
    prob = IntegralProblem(formated_integrand, lb, ub; nout=3, batch=Threads.nthreads())

    C = solve(prob, int_alg; reltol, abstol, maxiters)
    U = C[2] / C[1]
    θ^2 * ( C[3] / C[1] - U^2 )
end