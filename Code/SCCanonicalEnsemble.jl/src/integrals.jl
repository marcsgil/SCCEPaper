identity!(x,θ,par) = nothing
unit_jacobian_det(x,θ,par) = one(eltype(x))

function energy_output_integrals(sol, i, Xs, θ, par, H, H2)
    N = size(Xs, 1)

    Z = exp(sol[end, 1] - θ * H(view(Xs, :, i), par)) * √abs(extract_det_jac(sol.u[1], N))

    output = [Z, Z * H(sol[N+1:2N, 1], par)]

    regularize(output), false
end

function energy_integrals(θ, par, f!, H, lb, ub; 
    reltol=1e-2, abstol=1e-3, dif_eq_alg=nothing, int_alg=CubatureJLh(), 
    maxiters = 10^6, callback = caustic_callback, 
    change_of_variables! = identity!, jacobian_det = unit_jacobian_det)
    function formated_integrand(Xs, _)
        change_of_variables!(Xs,θ,par)
        jacobian_det(Xs,θ,par) * integrand(Xs, θ, par, f!, H, energy_output_integrals;
        dif_eq_alg, reltol, abstol, save_start=false, save_everystep=false, verbose=false, callback)
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
    maxiters = 10^6, callback = caustic_callback, 
    change_of_variables! = identity!, jacobian_det = unit_jacobian_det)
    function formated_integrand(Xs, _)
        change_of_variables!(Xs,θ,par)
        jacobian_det(Xs,θ,par) * integrand(Xs, θ, par, f!, H, heat_output_integrals;
        dif_eq_alg, reltol, abstol, save_start=false, save_everystep=false, verbose=false, callback)
    end

    prob = IntegralProblem(formated_integrand, lb, ub; nout=3, batch=Threads.nthreads())

    C = solve(prob, int_alg; reltol, abstol, maxiters)
    U = C[2] / C[1]
    θ^2 * ( C[3] / C[1] - U^2 )
end