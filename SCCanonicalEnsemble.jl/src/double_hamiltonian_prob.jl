function u0(X)
    #Given the initial center X, builds the initial conditions for the Double Hamilonian Problem
    T = eltype(X)
    N = length(X)
    vcat(zero(X), X, vec(vcat(zeros(eltype(X), N, N), I(N))), zero(T))
end

#=
The dimension of the phase space of a given vector u.
if length(X) == N, then phase_space_dim(u0(X)) = N.
=#
phase_space_dim(u) = Int((-1 + √(2length(u) - 1)) / 2)

#Exctracting the jacobian and its determinant of a vector u.
extract_jac_x(u, N) = view(reshape((@view u[2N+1:end-1]), 2N, N), N+1:2N, :)
extract_jac_x(u) = extract_jac_x(u, phase_space_dim(u))

extract_det_jac(u, N) = extract_jac_x(u, N) |> det
extract_det_jac(u) = extract_jac_x(u) |> det

#Callbacks to interrupt the trajectories upon crossing a caustic
function annul!(integrator)
    integrator.u = zero(integrator.u)
end

caustic_cross_contidion(u, t, integrator) = extract_det_jac(u) < 0 && extract_det_jac(integrator.uprev) > 0
caustic_callback = DiscreteCallback(caustic_cross_contidion, annul!, save_positions=(false, false))

function F!(du, u, par, f!, J, N)
    #This function defines the differential equation that needs to be solved

    #Differential equation for y and x
    f!(view(du, 1:2N), view(u, 1:2N), par)

    #Differential equation for the jacobians
    J.op.u .= view(u, 1:2N)

    for j in 1:N
        mul!(view(du, 2j*N+1:2*(j+1)*N), J, view(u, 2j*N+1:2*(j+1)*N))
    end

    #Differential equation for the area Δ
    du[end] = view(u, 1:N) ⋅ view(du, N+1:2N)

    nothing
end

function integrand(Xs, θ, par, f!, H, output_func; dif_eq_alg=nothing, kwargs...)
    u₀ = u0(view(Xs, :, 1))
    N = size(Xs, 1)

    J = JacVec((du, u) -> f!(du, u, par), view(u₀, 1:2N))

    prob = ODEProblem((du, u, par, t) -> F!(du, u, par, f!, J, N), u₀, (0, θ / 2), par)

    prob_func(prob, i, repeat) = remake(prob, u0=u0(view(Xs, :, i)))

    H2 = squared_hamiltonian_symbol(H, N ÷ 2)

    ensemble_prob = EnsembleProblem(prob; prob_func, output_func=(sol, i) -> output_func(sol, i, Xs, θ, par, H, H2))

    solve(ensemble_prob, dif_eq_alg, trajectories=size(Xs, 2); kwargs...) |> stack
end