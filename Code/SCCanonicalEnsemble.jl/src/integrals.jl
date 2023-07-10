function energy_output_integrals(sol, i, Xs, θ, par, H, H2)
    N = size(Xs, 1)

    Z = exp(sol[end, 1] - θ * H(view(Xs, :, i), par)) * √abs(extract_det_jac(sol.u[1], N))

    output = [Z, Z * H(sol[N+1:2N, 1], par)]

    isfinite(sum(output)) ? (output, false) : (zero(output), false)
end

function energy_integrals(θ, par, f!, H, ub, lb; reltol=1e-2, abstol=1e-3, dif_eq_alg=nothing, int_alg=CubatureJLh())
    formated_integrand(Xs, _) = integrand(Xs, θ, par, f!, H, energy_output_integrals;
        dif_eq_alg, reltol, abstol, save_start=false, save_everystep=false, verbose=false, callback=caustic_callback)

    prob = IntegralProblem(formated_integrand, lb, ub, nout=2, batch=Threads.nthreads())

    U = solve(prob, int_alg; reltol, abstol)
    U[2] / U[1]
end

function heat_output_integrals(sol, i, Xs, θ, par, H, H2)
    N = size(Xs, 1)

    Z = exp(sol[end, 1] - θ * H(view(Xs, :, i), par)) * √abs(extract_det_jac(sol.u[1], N))

    output = [Z, Z * H(sol[N+1:2N, 1], par), Z * H2(sol[N+1:2N, 1], par)]

    isfinite(sum(output)) ? (output, false) : (zero(output), false)
end

function heat_integrals(θ, par, f!, H, ub, lb; reltol=1e-2, abstol=1e-3, dif_eq_alg=nothing, int_alg=CubatureJLh())
    formated_integrand(Xs, _) = integrand(Xs, θ, par, f!, H, heat_output_integrals;
        dif_eq_alg, reltol, abstol, save_start=false, save_everystep=false, verbose=false, callback=caustic_callback)

    prob = IntegralProblem(formated_integrand, lb, ub, nout=3, batch=Threads.nthreads())

    C = solve(prob, int_alg; reltol, abstol)
    U = C[2] / C[1]
    U, θ^2 * ( C[3] / C[1] - U^2 )
end