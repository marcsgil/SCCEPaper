module ClassicalCanonicalEnsemble

export classical_partion_function,classical_energy,classical_heat

using Reexport
@reexport using Integrals,IntegralsCubature,IntegralsCuba
using ThreadsX

transform(t) = t/(1-t^2)
g(t) = (1+t^2)/(1-t^2)^2

function classical_partion_function(θ::Real,H,d::Integer,par=nothing)

    function integrand!(y,u,par)
        Threads.@threads for i in axes(u,2)
            y[i] = exp(-θ*H(transform.(view(u,:,i)),par))*prod(g,view(u,:,i))
        end
    end

    prob = IntegralProblem(integrand!,fill(-1,2d),fill(1,2d),par,batch=2)
    solve(prob,CubatureJLh()).u[1]/(2π)^d
end

function classical_partion_function(θs::AbstractArray,H,d::Integer,par=nothing; indicator_function = (x,par) -> true)
    ThreadsX.map(θ->classical_partion_function(θ,H,d,par, indicator_function = indicator_function),θs)
end

function classical_energy(θ,d::Integer,par=nothing; hamiltonian=nothing,potential=nothing, alg = CubatureJLh(),integration_limit=50,kwargs...)

    @assert !isnothing(hamiltonian) || !isnothing(potential) "Please inform the hamiltonian or the potential through kwargs."

    @assert isnothing(hamiltonian) || isnothing(potential) "Please inform only the hamiltonian or only the potential."

    f = isnothing(hamiltonian) ? potential : hamiltonian

    function integrand!(y,u,par)
        Threads.@threads for i in axes(u,2)
            E = f(view(u,:,i),par)
            y[1,i] = exp(-θ*E)
            y[2,i] = E*y[1,i]
        end
    end

    prob = IntegralProblem(integrand!,fill(-integration_limit,2d),fill(integration_limit,2d),par,batch=2,nout=2)
    sol = solve(prob,alg;kwargs...).u
    isnothing(hamiltonian) ? sol[2]/sol[1] + d/(2θ) : sol[2]/sol[1]
end

function classical_energy(θs::AbstractArray,d::Integer,par=nothing;  kwargs...)
    ThreadsX.map(θ->classical_energy(θ,d,par; kwargs...),θs)
end

function classical_heat(θ,d::Integer,par=nothing; hamiltonian=nothing,potential=nothing,integration_limit=5,alg=CubatureJLh(), kwargs...)

    @assert !isnothing(hamiltonian) || !isnothing(potential) "Please inform the hamiltonian or the potential through kwargs."

    @assert isnothing(hamiltonian) || isnothing(potential) "Please inform only the hamiltonian or only the potential."

    f = isnothing(hamiltonian) ? potential : hamiltonian

    function integrand!(y,u,par)
        Threads.@threads for i in axes(u,2)
            E = f(view(u,:,i),par)
            y[1,i] = exp(-θ*E)
            y[2,i] = E*y[1,i]
            y[3,i] = E*y[2,i]
        end
    end

    prob = IntegralProblem(integrand!,fill(-integration_limit,2d),fill(integration_limit,2d),par,batch=2,nout=3)
    sol = solve(prob,alg;kwargs...).u
    isnothing(hamiltonian) ? θ^2 * ( sol[3]/sol[1] - ( sol[2]/sol[1] )^2 ) + d/2 : θ^2 * ( sol[3]/sol[1] - ( sol[2]/sol[1] )^2 )
end

function classical_heat(θs::AbstractArray,d::Integer,par=nothing; kwargs...)
    ThreadsX.map(θ->classical_heat(θ,d,par;kwargs...),θs)
end

end