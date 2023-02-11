module ClassicalCanonicalEnsemble

export classical_partion_function,classical_energy

using Integrals,IntegralsCubature,ForwardDiff,ThreadsX

transform(t) = t/(1-t^2)
g(t) = (1+t^2)/(1-t^2)^2

function classical_partion_function(θ,H,d::Integer,par=nothing; indicator_function = (x,par) -> true)

    function integrand(x,par)
        X = transform.(x)
        if indicator_function(X,par)
            exp(-θ*H(X,par))*prod(g,x)
        else
            zero(eltype(x))
        end
    end

    function integrand!(y,u,par)
        Threads.@threads for i in axes(u,2)
            y[i] = integrand(view(u,:,i),par)
        end
    end

    prob = IntegralProblem(integrand!,fill(-1,2d),fill(1,2d),par,batch=2)
    sol = solve(prob,CubatureJLh()).u[1]/(2π)^d
end

function classical_energy(θ,H,d::Integer,par=nothing; indicator_function = (x,par) -> true)
    function integrand(x,par)
        X = transform.(x)
        if indicator_function(X,par)
            exp(-θ*H(X,par))*prod(g,x)
        else
            zero(eltype(x))
        end
    end

    function integrand!(y,u,par)
        Threads.@threads for i in axes(u,2)
            y[1,i] = integrand(view(u,:,i),par)
            y[2,i] = integrand(view(u,:,i),par)*H(transform.(view(u,:,i)),par)
        end
    end

    prob = IntegralProblem(integrand!,fill(-1,2d),fill(1,2d),par,batch=2,nout=2)
    sol = solve(prob,CubatureJLh()).u
    sol[2]/sol[1]
end


end