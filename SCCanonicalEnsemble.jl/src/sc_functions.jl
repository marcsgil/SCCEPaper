function Z_integrand(u,energy_term)
    exp( u[end] - energy_term )*√abs(extract_det_jac(u))
end

function energy_output(sol,i,θs,par,nodes,weights,H)
    N=phase_space_dim(sol[end])
    output = Array{eltype(θs)}(undef,2,length(θs))

    E = H(nodes[i],par)
    for (n,u) in enumerate(sol.u)
        output[1,n] = weights[i]*Z_integrand(u,θs[n]*E)
        output[2,n] = H(view(u,N+1:2N),par)
    end
    output,false
end

function Z_integrandMonteCarlo(u,θ,par,nodes,i)
    exp(u[end])*√abs(extract_det_jac(u))
end

function energy_outputMonteCarlo(sol,i,θ,par,nodes,weights,H)
    N=phase_space_dim(sol.u[end])
    [Z_integrandMonteCarlo(sol.u[end],θ,par,nodes[i],i),H(view(sol.u[end],N+1:2N),par)],false
end

function energy_reduction(sols,θ)
    if ndims(sols) == 2
        sum(prod,eachcol(sols))/sum(first,view(sols,1,:))
    elseif ndims(sols) == 3
        [ energy_reduction(sol,θ) for sol in eachslice(sols,dims=2)  ]
    else
        @error "Unsuported number of dimensions."
        sols
    end
end

function sample_points(θ,par,N,H,d)
    xs = metropolisHastings(x -> exp(-θ*H(x,par) ),N,2d)
    [xs[:,n] for n in 1:N], nothing
end

function energyMonteCarlo(θ,par,H,d,f,N;alg=BS3(),reltol=1e-1,abstol=1e-2,callback=strong_callback)
    solve_equations(θ,par,f,(θ,par)->sample_points(θ,par,N,H,d),H,
    output_func=energy_outputMonteCarlo,reduction=energy_reduction,alg=alg,reltol=reltol,abstol=abstol,callback=callback)
end

function analysis_output(sol,i,θs,par,nodes,weights,H)
    output = Array{eltype(θs)}(undef,3,length(θs))
    E = H(nodes[i],par)
    for (n,u) in enumerate(sol.u)
        output[2,n] = extract_det_jac(u) 
        output[1,n] = exp( u[end] - θs[n]*E )*√abs(output[2,n])
        output[3,n] = u[end]
    end
    output,false
end

function analysis_outputMC(sol,i,θs,par,nodes,weights)
    output = Array{eltype(θs)}(undef,3,length(θs))
    E = H(nodes[i],par)
    for (n,u) in enumerate(sol.u)
        output[2,n] = extract_det_jac(u) 
        output[1,n] = exp( u[end] )*√abs(output[2,n])
        output[3,n] = u[end]
    end
    output,false
end