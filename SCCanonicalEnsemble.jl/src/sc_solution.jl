function u0(X)
    #Given the initial center X, builds the initial conditions for the initial value problem.
    T = eltype(X)
    N = length(X)
    vcat(zero(X),X,vec(vcat(zeros(eltype(X),N,N),I(N))),zero(T))
end

function extract_jac_x(u)
    N = Int((-1+√(2length(u)-1))/2)
    view(reshape((@view u[2N+1:end-1]),2N,N),N+1:2N,:)
end

function extract_det_jac(u)
    N = Int((-1+√(2length(u)-1))/2)
    view(reshape((@view u[2N+1:end-1]),2N,N),N+1:2N,:) |> det
end

function annul!(integrator)
    integrator.u = zero(integrator.u)
end

area_condition(u,t,integrator) = u[end] > 0 && integrator.uprev[end] < 0
caustic_cross_contidion(u,t,integrator) = extract_det_jac(u)
caustic_callback = ContinuousCallback(caustic_cross_contidion,annul!,save_positions=(false,false))
disc_caustic_cross_contidion(u,t,integrator) = extract_det_jac(u) < 0 && extract_det_jac(integrator.uprev) > 0
disc_caustic_callback = DiscreteCallback(disc_caustic_cross_contidion,annul!,save_positions=(false,false))
area_callback = DiscreteCallback(area_condition,annul!,save_positions=(false,false))
strong_callback = CallbackSet(disc_caustic_callback,area_callback)

function phase_space_dim(u)
    Int((-1+√(2length(u)-1))/2)
end

function solve_equations(θ,par,f!,getNodesAndWeights,H;
    output_func=(sol,i,θ,par,node,weight,H)->(sol,false),reduction=(sols,θ)->sols,alg=BS3(),
    reltol=1e-1,abstol=1e-2,callback=caustic_callback)
    #=Returns the solution of the ModifiedHamiltonianProblem 
    for each initial condition in nodes, in the interval (0,θ_max)=#

    θ = float.(θ)
    
    nodes,weights = applicable(getNodesAndWeights,par) ? getNodesAndWeights(par) : getNodesAndWeights(θ,par)

    initial_condition = u0(nodes[1])
    N = phase_space_dim(initial_condition)
    J = JacVec((du,u) -> f!(du,u,par),view(initial_condition,1:2N))

    function F!(du,u,par,J,N)
        #Differential equation for y and x
        f!(view(du,1:2N),view(u,1:2N),par)
    
        #Differential equation for the jacobians
        J.op.u .= view(u,1:2N)
    
        for j in 1:N
            mul!(view(du, 2j*N+1:2*(j+1)*N ),J,view(u,2j*N+1:2*(j+1)*N ))
        end
    
        #Differential equation for the area Δ
        du[end] = view(u,1:N) ⋅ view(du,N+1:2N)
    
        nothing
    end

    prob = ODEProblem((du,u,par,t)->F!(du,u,par,J,N),initial_condition,(0,last(θ)/2),par)

    #Changes the initial condition after each run
    function prob_func(prob,i,repeat)
        remake(prob,u0=u0(nodes[i]))
    end

    H2 = squared_hamiltonian_symbol(H,N ÷ 2)

    formated_ouput_func(sol,i) = applicable(output_func,sol,1,first(θ),par,nodes,weights,H) ? 
    output_func(sol,i,θ,par,nodes,weights,H) : output_func(sol,i,θ,par,nodes,weights,H,H2)    

    ensemble_prob = EnsembleProblem(prob,prob_func=prob_func,output_func=formated_ouput_func)

    if θ isa AbstractArray
        sols = solve(ensemble_prob,alg,trajectories=length(nodes),reltol=reltol,abstol=abstol,
        callback=callback,saveat=θ/2)
    else
        sols = solve(ensemble_prob,alg,trajectories=length(nodes),reltol=reltol,abstol=abstol,
        callback=callback,save_start=false,save_everystep=false)
    end   

    reduction(sols,θ)
end