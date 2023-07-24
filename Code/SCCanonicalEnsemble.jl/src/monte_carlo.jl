function generate_mc_samples(N, θ, par, H, d)
    density(x) = -θ * H(x, par)
    model = DensityModel(density)
    proposal = MetropolisHastings(RandomWalkProposal(MvNormal(zeros(2d), I)))
    samples = sample(model, proposal, N; chain_type=Chains)
    view(samples.value.data, :, 1:2d, 1) |> permutedims
end

function energy_output_mc(sol, i, Xs, θ, par, H, H2)
    N = size(Xs, 1)

    Z = exp(sol[end, 1]) * √abs(extract_det_jac(sol.u[1], N))
    output = [Z, Z*H]

    regularize(output), false
end

function energy_mc(θ, par, f!, H, d;
    reltol=1e-2, abstol=1e-3, dif_eq_alg=nothing, batchsize=2^12,
    maxiters=10^5, callback=caustic_callback)

    @assert batchsize ≤ maxiters

    s = Series(Mean(), Variance())

    @showprogress for _ ∈ 1:maxiters÷batchsize
        Xs = generate_mc_samples(batchsize, θ, par, H, d)
        estimators = integrand(Xs, θ, par, f!, H, energy_output_mc;
            dif_eq_alg, reltol, abstol, save_start=false, save_everystep=false, verbose=false, callback)

        Z_estimator = view(estimators, 1, :) |> mean 
        U_estimator = view(estimators, 2, :) |> mean

        OnlineStats.fit!(s, U_estimator / Z_estimator )

        val, σ2 = value.(s.stats)
        if abs(σ2 / val) < reltol || √σ2 < abstol
            #break
        end
    end

    val, σ2 = value.(s.stats)
    val, √σ2
end

function heat_output_mc(sol, i, Xs, θ, par, H, H2)
    N = size(Xs, 1)

    Z = exp(sol[end, 1]) * √abs(extract_det_jac(sol.u[1], N))
    output = [Z, Z * H(sol[N+1:2N, 1], par), Z * H2(sol[N+1:2N, 1], par)]

    regularize(output), false
end

function heat_mc(θ, par, f!, H, d;
    reltol=1e-2, abstol=1e-3, dif_eq_alg=nothing, batchsize=2^12,
    maxiters=10^5, callback=caustic_callback)

    @assert batchsize ≤ maxiters

    s = Series(Mean(), Variance())

    @showprogress for _ ∈ 1:maxiters÷batchsize
        Xs = generate_mc_samples(batchsize, θ, par, H, d)
        estimators = integrand(Xs, θ, par, f!, H, heat_output_mc;
            dif_eq_alg, reltol, abstol, save_start=false, save_everystep=false, verbose=false, callback)

        inv_Z_estimator = view(estimators, 1, :) |> mean |> inv
        U_estimator = view(estimators, 2, :) |> mean
        U2_estimator = view(estimators, 3, :) |> mean
        c_estimator = θ^2 * (U2_estimator - U_estimator^2*inv_Z_estimator)*inv_Z_estimator

        OnlineStats.fit!(s, c_estimator)

        val, σ2 = value.(s.stats)
        if abs(σ2 / val) < reltol || √σ2 < abstol
            #break
        end
    end

    val, σ2 = value.(s.stats)
    val
end
