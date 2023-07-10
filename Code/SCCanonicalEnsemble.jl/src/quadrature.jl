function energy_output_quadrature(sol, i, Xs, θs, par, H, H2)
    N = size(Xs,1)
    output = Array{eltype(θs)}(undef,2,length(θs))

    for (n,u) in enumerate(sol.u)
        output[1,n] = exp(u[end] - θs[n] * H(view(Xs, :, i), par)) * √abs(extract_det_jac(u, N))
        output[2,n] = output[1,n] * H(view(u,N+1:2N),par)
    end
    isfinite(sum(output)) ? (output, false) : (zero(output), false)
end

function energy_quadrature(θ::AbstractArray, par, f!, H, Xs, ws; reltol=1e-2, abstol=1e-3, dif_eq_alg=nothing)
    integrands = integrand(Xs, θ, par, f!, H, energy_output_quadrature;
        dif_eq_alg, reltol, abstol, saveat=θ/2, verbose=false, callback=caustic_callback)
    
    U = mapreduce( (x,w) -> w .* x,+, eachslice(integrands,dims=3), ws )
    @views U[2,:] ./ U[1,:]
end