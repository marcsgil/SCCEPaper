function wigner_output(sol, i, Xs, θ, par, H, H2)
    N = size(Xs, 1)

    Z = exp(sol[end, 1] - θ * H(view(Xs, :, i), par)) * √abs(extract_det_jac(sol.u[1], N))

    regularize(Z), false
end