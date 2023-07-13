using DynamicalSystems, SCCanonicalEnsemble

function nelson!(du,u,μ,t)
    du[1] = -( 2μ*u[3] + 2 * (u[3]^2/2-u[4]) * u[3] )
    du[2] = - ((u[4]-u[3]^2/2))
    du[3] = u[1]
    du[4] = u[2]
end

V(q, μ) = (q[1]^2 / 2 - q[2])^2 + μ * q[1]^2
H(x, μ) = (x[1]^2 + x[2]^2) / 2 + (x[3]^2 / 2 - x[4])^2 + μ * x[3]^2
f! = get_equations_of_motion(H, 2)
g!(du,u,par,t) = f!(du,u,par)
##
x0 =  rand(4)
u0 = vcat(zeros(4), x0)
tspan = (0.0, 10.0)
μ = 1

prob = ODEProblem(g!, u0, tspan, μ)
system = ContinuousDynamicalSystem(prob)
lyapunovspectrum(system, 10^4)