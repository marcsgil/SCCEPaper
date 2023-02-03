#=function J(d)
    J = zeros(Int,2d,2d)
    for j in axes(J,2)
        for i in axes(J,1)
            if i+d == j
                J[i,j] = -1
            elseif j+d == i
                J[i,j] = 1
            end
        end
    end
    J
end=#

function get_equations_of_motion2(H,d)
    Meta.parse( "@variables " * prod("u$n::Real " for n in 1:4d) * "par::Real t::Real" ) |> eval

    u = Meta.parse("[" * prod(["u$n, " for n in 1:4d]) * "]") |> eval
    y = Meta.parse("[" * prod(["u$n, " for n in 1:2d]) * "]") |> eval
    x = Meta.parse("[" * prod(["u$n, " for n in 2d+1:4d]) * "]") |> eval

    ℍ = H(x+im*.5*J(d)*y,par) + H(x-im*.5*J(d)*y,par) |> real

    simplify.(build_function( vcat(Symbolics.gradient(-ℍ, x),Symbolics.gradient( ℍ, y),sum(y .* Symbolics.gradient( ℍ, y))), u,par,t))[1] |> eval
end

#=area_condition(u,t,integrator) = u.Δ > 0 && integrator.uprev.Δ < 0
caustic_cross_contidion(u,t,integrator) = det(u.jac_x)
caustic_callback = ContinuousCallback(caustic_cross_contidion,affect!,save_positions=(false,false))
area_callback = DiscreteCallback(area_condition,affect!,save_positions=(false,false))
strong_callback = CallbackSet(caustic_callback,area_callback)=#

u02(X) = vcat(vcat(zero(X),X),zero(eltype(X)))


function solve_equations2(θ,par,f,getNodesAndWeights,H;
    output_func=(sol,i,θ,par,node,weight,H)->(sol,false),reduction=(sols,θ)->sols,alg=BS3(),
    reltol=1e-1,abstol=1e-2,callback=strong_callback)
    #=Returns the solution of the ModifiedHamiltonianProblem 
    for each initial condition in nodes, in the interval (0,θ_max)=#

    θ = float.(θ)
    
    nodes,weights = applicable(getNodesAndWeights,par) ? getNodesAndWeights(par) : getNodesAndWeights(θ,par)
    prob = ODEProblem(f,u02(nodes[1]),(0,last(θ)/2),par)

    #Changes the initial condition after each run
    function prob_func(prob,i,repeat)
        remake(prob,u0=u02(nodes[i]))
    end

    ensemble_prob = EnsembleProblem(prob,prob_func=prob_func,output_func=(sol,i)->output_func(sol,i,θ,par,nodes,weights,H))

    sols = solve(ensemble_prob,alg,trajectories=length(nodes),reltol=reltol,abstol=abstol,
    callback=callback,save_start=false,save_everystep=false)

    reduction(sols,θ)
end