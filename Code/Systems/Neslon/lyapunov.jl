using DynamicalSystems, SCCanonicalEnsemble
include("../../plot_config.jl")

function nelson!(du, u, μ, t)
    du[1] = -(2μ * u[3] + 2 * (u[3]^2 / 2 - u[4]) * u[3])
    du[2] = -((u[4] - u[3]^2 / 2))
    du[3] = u[1]
    du[4] = u[2]
end

V(q, μ) = (q[1]^2 / 2 - q[2])^2 + μ * q[1]^2
H(x, μ) = (x[1]^2 + x[2]^2) / 2 + (x[3]^2 / 2 - x[4])^2 + μ * x[3]^2
f! = get_equations_of_motion(H, 2)
g!(du, u, par, t) = f!(du, u, par)
##
x0 = rand(4)
u0 = vcat(zeros(4), x0)
tspan = (0.0, 10.0)
μ = 0.6

prob = ODEProblem(nelson!, x0, tspan, μ)
system = ContinuousDynamicalSystem(prob)
lyapunovspectrum(system, 10^4)
##
x1 = [1.1, 1.9, 1.8, 0.2]
x2 = [4.157092673031904, 4.289582071149792, 3.618293817815962, 3.7340834096292492]
x3 = [0.2919906617580341, 0.7477332931666828, 0.7120035011399841, 0.5393250435426595,]
##
x0 = rand(4)
μs = LinRange(0.05, 2, 100)
λs = similar(μs)

Threads.@threads for n ∈ eachindex(μs)
    tspan = (0.0, 10.0)
    prob = ODEProblem(nelson!, x0, tspan, μs[n])
    system = ContinuousDynamicalSystem(prob)
    λs[n] = lyapunov(system, 10^6)
end

λs

lines(μs, λs)