using GLMakie, DynamicalSystems, StaticArrays, ForwardDiff
using OrdinaryDiffEq: Vern9

diffeq = (alg=Vern9(), abstol=1e-9, reltol=1e-9)

function henonheiles_rule(u, p, t)
    du1 = u[3]
    du2 = u[4]
    du3 = -u[1] - 2.0 * u[1] * u[2]
    du4 = -u[2] - (u[1]^2 - u[2]^2)
    return SVector(du1, du2, du3, du4)
end

V(x, y, μ) = μ * x^2 + (x^2 / 2 - y)^2
H(x, y, px, py, μ) = (px^2 + py^2) / 2 + V(x, y, μ)

function nelson(u, μ, t)
    du1 = ForwardDiff.derivative(px -> H(u[1], u[2], px, u[4], μ), u[3])
    du2 = ForwardDiff.derivative(py -> H(u[1], u[2], u[3], py, μ), u[4])
    du3 = -ForwardDiff.derivative(x -> H(x, u[2], u[3], u[4], μ), u[1])
    du4 = -ForwardDiff.derivative(y -> H(u[1], y, u[3], u[4], μ), u[2])
    return SVector(du1, du2, du3, du4)
end
##
#H(u,μ) = (u[3]^2 + u[4]^2) / 2 + μ * u[1]^2 + (u[1]^2/2 - u[2])^2
μ = 2.0
ds = CoupledODEs(nelson, [0.0, -0.25, 0.42081, 0.0], μ; diffeq)
##
E = 4.8
u0s = Vector{Float64}[]

for _ ∈ 1:20
    for _ ∈ 1:1000
        x = randn()
        y = randn()
        K = E - V(x, y, μ)
        K>0 ? break : nothing
    end
    px = (rand() - 0.5) * 2√(2K)
    py = √(2K - px^2)
    push!(u0s,[x,y,px,py])
end
##
trs = [trajectory(ds, 10000, u0)[1][:, SVector(1, 2, 3)] for u0 ∈ u0s]
j = 2 # the dimension of the plane

figure, ax3D, ax2D = interactive_poincaresos_scan(trs, j;
    linekw=(transparency=true,), scatterkw=(alpha=.8, markersize=4,))
##
u0 = 10 * randn(4)
tr = trajectory(ds, 10000, u0)[1]
map(n -> H(tr[n], 2), eachindex(tr))
##

##
ds = PredefinedDynamicalSystems.henonheiles()
ds = CoupledODEs(ds, diffeq)

u0s = [
    [0.0, -0.25, 0.42081, 0.0],
    [0.0, 0.1, 0.5, 0.0],
    [0.0, -0.31596, 0.354461, 0.0591255]
]

#u0s = [randn(4) for _ ∈ 1:3]
# inputs
trs = [trajectory(ds, 10000, u0)[1][:, SVector(1, 2, 3)] for u0 ∈ u0s]
j = 2 # the dimension of the plane

interactive_poincaresos_scan(trs, j; linekw=(transparency=true,))
##
propertynames(ds)
getproperty(ds, :diffeq)
PredefinedDynamicalSystems.henonheiles_rule