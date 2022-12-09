using DifferentialEquations,ComponentArrays,SparseDiffTools,LinearAlgebra
using Parameters: @unpack

function u0(X)
    #Given the initial center X, builds the initial conditions for the initial value problem.
    temp = X*transpose(X)
    ComponentArray( y=zero(X),x=X,Δ=zero(eltype(X)),jac_y= zero(temp),jac_x=one(temp))
end

function ModifiedHamiltonianProblem(fy,fx,initial_condition,tspan,par=nothing)
    function F!(du,u,par,t)
        @unpack x,y,Δ,jac_y,jac_x = u
    
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

function calculate_expectation(θ,fy,fx,(nodes,weights),output_func,reduction,par=nothing)
    #=Returns the solution of the ModifiedHamiltonianProblem 
    for each initial condition in nodes, in the interval (0,θ_max)=#
    
    prob = ModifiedHamiltonianProblem(fy,fx,u0(nodes[1]),(0,θ/2),par)

    #Changes the initial condition after each run
    function prob_func(prob,i,repeat)
        remake(prob,u0=ComponentArray(prob.u0,x=nodes[i]))
    end

    ensemble_prob = EnsembleProblem(prob,prob_func=prob_func,output_func=(sol,i)->energy_output(sol,i,θ,par,(nodes,weights),H))

    condition(u,t,integrator) = det(u.jac_x)           #This will return true if det(u.jac_x)==0
    affect!(integrator) = terminate!(integrator)       #We terminate the solution when the condition is true
    cb = ContinuousCallback(condition,affect!,save_positions=(true,false))

    sols = solve( ensemble_prob,BS3(),trajectories=length(nodes),reltol=1e-1,abstol=1e-2,save_everystep=false,save_start=false,callback=cb ) 

    reduction(sols)
end