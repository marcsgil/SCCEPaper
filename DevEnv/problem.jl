using DifferentialEquations,SparseDiffTools

function F!(du,u,par,J,N)
    #Differential equation for y and x
    g!(view(du,1:2N),view(u,1:2N),par)

    #Differential equation for the jacobians
    J.op.u .= view(u,1:2N)

    for j in 1:N
        mul!(view(du, 2j*N+1:2*(j+1)*N ),J,view(u,2j*N+1:2*(j+1)*N ))
    end

    #Differential equation for the area Δ
    du[end] = view(u,1:N) ⋅ view(du,N+1:2N)

    nothing
end

struct DoubleHamiltonianProblem{uType, tType}
    prob::ODEProblem
    θ::tType
    J::FunctionOperator

    function MarkovProblem(f!, X, tspan, p, dt; kwargs...)
        N = length(X)
        T = eltype(X)
        u₀ = vcat(zero(X),X,vec(vcat(zeros(eltype(X),N,N),I(N))),zero(T))

        J = JacVec((du,u) -> f!(du,u,par),view(u₀,1:2N))

        prob = ODEProblem((du,u,p,t) -> f(du,u,p,t,dt), u0, (0,θ/2), p; kwargs...)
        new{T, float(typeof(θ))}(prob, dt, rng)
    end
end

struct MarkovSolution{uType, tType}
    sol::ODESolution
    function MarkovSolution(sol)
        new{eltype(eltype(sol.u)), eltype(sol.t)}(sol)
    end
end

function solve(prob::MarkovProblem; kwargs...)
    sol = solve(prob.prob, FunctionMap(), dt=prob.dt, kwargs...)
    MarkovSolution(sol)
end

function Base.show(io::IO, sol::MarkovSolution)
    A = sol.sol
    println(io, string("t: ", A.t))
    print(io, "u: ")
    show(io, A.u)
end

function Base.show(io::IO, m::MIME"text/plain", sol::MarkovSolution)
    A = sol.sol
    println(io, string("t: ", A.t))
    print(io, "u: ")
    show(io, m, A.u)
end

function plot(sol::MarkovSolution; kwargs...)
    plot(sol.sol; kwargs...)
end

## Test

using Random
using Distributions
using Plots

function sir_markov!(du,u,p,t,dt,rng)
    (S,I) = u
    (β,γ) = p
    ifrac = 1-exp(-β*I*dt)
    rfrac = 1-exp(-γ*dt)
    infection=rand(rng, Binomial(S,ifrac))
    recovery=rand(rng, Binomial(I,rfrac))
    du[1] = S-infection
    du[2] = I+infection-recovery
    nothing
end

tspan = (0.0, 40.0)
dt = 0.1
u0 = [990, 10] # S,I
p = [0.0005, 0.25] # β,γ
rng = Xoshiro(1234)

prob = MarkovProblem(sir_markov!, u0, tspan, p, dt, rng)
sol = solve(prob)
plot(sol)