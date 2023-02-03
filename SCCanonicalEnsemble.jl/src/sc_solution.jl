function u0(X)
    #Given the initial center X, builds the initial conditions for the initial value problem.
    temp = X*transpose(X)
    ComponentArray( y=zero(X),x=X,Δ=zero(eltype(X)),jac_y= zero(temp),jac_x=one(temp))
end

function ModifiedHamiltonianProblem(fy,fx,initial_condition,tspan,par=nothing)
    function F!(du,u,par,t)
        Parameters.@unpack x,y,Δ,jac_y,jac_x = u
    
        #Differential equation for y and x
        du.y = fy(y,x,par)
        du.x = fx(y,x,par)

        #Differential equation for the area Δ
        du.Δ = y⋅du.x
    
        #Differential equation for the jacobians
        for n in axes(jac_y,2)
            du.jac_y[:,n] = @views auto_jacvec(y->fy(y,x,par), y, jac_y[:,n]) + auto_jacvec(x->fy(y,x,par), x, jac_x[:,n])
            du.jac_x[:,n] = @views auto_jacvec(y->fx(y,x,par), y, jac_y[:,n]) + auto_jacvec(x->fx(y,x,par), x, jac_x[:,n])
        end
        
    end
    
    ODEProblem( F!, initial_condition, tspan, par)
end

function affect!(integrator)
    integrator.u = zero(integrator.u)
    terminate!(integrator)
end

area_condition(u,t,integrator) = u.Δ > 0 && integrator.uprev.Δ < 0
caustic_cross_contidion(u,t,integrator) = det(u.jac_x)
caustic_callback = ContinuousCallback(caustic_cross_contidion,affect!,save_positions=(false,false))
area_callback = DiscreteCallback(area_condition,affect!,save_positions=(false,false))
strong_callback = CallbackSet(caustic_callback,area_callback)


function solve_equations(θ,par,fy,fx,getNodesAndWeights,H;
    output_func=(sol,i,θ,par,node,weight,H)->(sol,false),reduction=(sols,θ)->sols,alg=BS3(),
    reltol=1e-1,abstol=1e-2,callback=strong_callback)
    #=Returns the solution of the ModifiedHamiltonianProblem 
    for each initial condition in nodes, in the interval (0,θ_max)=#

    θ = float.(θ)
    
    nodes,weights = applicable(getNodesAndWeights,par) ? getNodesAndWeights(par) : getNodesAndWeights(θ,par)
    prob = ModifiedHamiltonianProblem(fy,fx,u0(nodes[1]),(0,last(θ)/2),par)

    #Changes the initial condition after each run
    function prob_func(prob,i,repeat)
        remake(prob,u0=ComponentArray(prob.u0,x=nodes[i]))
    end

    ensemble_prob = EnsembleProblem(prob,prob_func=prob_func,output_func=(sol,i)->output_func(sol,i,θ,par,nodes,weights,H))

    if typeof(θ) <: AbstractArray
        sols = solve(ensemble_prob,alg,trajectories=length(nodes),reltol=reltol,abstol=abstol,
        callback=callback,saveat=θ/2)
    else
        sols = solve(ensemble_prob,alg,trajectories=length(nodes),reltol=reltol,abstol=abstol,
        callback=callback,save_start=false,save_everystep=false)
    end   

    reduction(sols,θ)
end