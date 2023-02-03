function J(d)
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
end

function get_equations_of_motion(H,d)
    Meta.parse( "@variables " * prod("x$n::Real " for n in 1:2d) * prod("y$n::Real " for n in 1:2d)*"par::Real" ) |> eval

    exp_x = Meta.parse("[" * prod(["x$n, " for n in 1:2d]) * "]")
    exp_y = Meta.parse("[" * prod(["y$n, " for n in 1:2d]) * "]")

    x = eval(exp_x)
    y = eval(exp_y)

    ℍ = H(x+im*J(d)*y/2,par) + H(x-im*J(d)*y/2,par) |> real
    build_function( Symbolics.gradient(-ℍ, eval(exp_x)), y,x,par)[1] |> eval, build_function( Symbolics.gradient(ℍ, eval(exp_y)), y,x,par)[1] |> eval
end

function Z_integrand(u,energy_term)
    exp( u.Δ - energy_term )*√abs(det(u.jac_x))
end

function energy_output(sol,i,θs,par,nodes,weights,H)
    output = Array{eltype(θs)}(undef,2,length(θs))

    E = H(nodes[i],par)
    for (n,u) in enumerate(sol.u)
        output[1,n] = weights[i]*Z_integrand(u,θs[n]*E)
        output[2,n] = H(u.x,par)
    end
    output,false
end

function Z_integrandMonteCarlo(u,θ,par,nodes,i)
    exp(u.Δ)*√abs(det(u.jac_x))
end

function energy_outputMonteCarlo(sol,i,θ,par,nodes,weights,H)
    [Z_integrandMonteCarlo(sol[end],θ,par,nodes[i],i),H(sol[end].x,par)],false
end

function energy_reduction(sols,θ)
    if ndims(sols) == 2
        sum(prod,eachcol(sols))/sum(first,view(sols,1,:))
    elseif ndims(sols) == 3
        [ energy_reduction(sol,θ) for sol in eachslice(sols,dims=2)  ]
    else
        @error "Unsuported number of dimensions."
    end
end

function sample_points(θ,par,N,H,d)
    xs = metropolisHastings(x -> exp(-θ*H(x,par) ),N,2d)
    [SVector{2d}(xs[:,n]) for n in 1:N], nothing
end

function energyMonteCarlo(θ,par,H,d,fy,fx,N;alg=BS3(),reltol=1e-1,abstol=1e-2,callback=strong_callback)
    solve_equations(θ,par,fy,fx,(θ,par)->sample_points(θ,par,N,H,d),H,
    output_func=energy_outputMonteCarlo,reduction=energy_reduction,alg=alg,reltol=reltol,abstol=abstol,callback=callback)
end

function analysis_output(sol,i,θs,par,nodes,weights,H)
    output = Array{eltype(θs)}(undef,3,length(θs))
    E = H(nodes[i],par)
    for (n,u) in enumerate(sol.u)
        output[1,n] = exp( u.Δ - θs[n]*E )*√abs(det(u.jac_x))
        output[2,n] = det(u.jac_x)
        output[3,n] = u.Δ
    end
    output,false
end

function analysis_outputMC(sol,i,θs,par,nodes,weights)
    output = Array{eltype(θs)}(undef,3,length(θs))
    E = H(nodes[i],par)
    for (n,u) in enumerate(sol.u)
        output[1,n] = exp( u.Δ )*√abs(det(u.jac_x))
        output[2,n] = det(u.jac_x)
        output[3,n] = u.Δ
    end
    output,false
end