function integrand(Xs,θ,par,f!,H,output_func)
    u₀ = u0(view(Xs,:,1))
    N = size(Xs,1)

    J = JacVec((du,u) -> f!(du,u,par),view(u₀,1:2N))

    prob = ODEProblem((du,u,par,t)->F!(du,u,par,f!,J,N),u₀,(0,θ/2),par)

    function prob_func(prob,i,repeat)
        remake(prob,u0=u0(view(Xs,:,i)))
    end

    ensemble_prob = EnsembleProblem(prob;prob_func,output_func=(sol,i)->output_func(sol,i,Xs,H,θ,par,N))

    solve(ensemble_prob,BS3(),trajectories=size(Xs,2),reltol=1e-2,abstol=1e-3,
    save_everystep=false,save_start=false,callback=disc_caustic_callback,verbose=false) |> stack
end

function energy_output_cubature(sol,i,Xs,H,θ,par,N)
    Z = exp( sol[end,1] - θ*H(view(Xs,:,i),par) )*√abs(extract_det_jac(sol.u[1])) 
    output = [Z,Z*H(sol[N+1:2N,1],par)]
    isfinite(sum(output)) ? (output,false) : (zero(output),false)
end

function energyCubature(θ,par,f!,H,ndims)
    ub = [10. for n ∈ 1:ndims]
    lb = -ub
    formated_integrand(Xs,_) = integrand(Xs,θ,par,f!,H,energy_output_cubature)
    prob = IntegralProblem(formated_integrand, lb, ub, nout = 2, batch = 64)
    U = solve(prob, CubatureJLh(), reltol = 1e-2, abstol = 1e-3)
    U[2] / U[1]
end