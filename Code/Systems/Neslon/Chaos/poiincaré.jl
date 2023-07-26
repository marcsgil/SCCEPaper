using GLMakie, DynamicalSystems, StaticArrays
using OrdinaryDiffEq: Vern9, ODEProblem

diffeq = (alg = Vern9(), abstol = 1e-9, reltol = 1e-9)

function nelson(u,μ,t)
    du1 = -( 2μ*u[3] + 2 * (u[3]^2/2-u[4]) * u[3] )
    du2 = - ((u[4]-u[3]^2/2))
    du3 = u[1]
    du4 = u[2]
    @SVector [du1,du2,du3,du4]
end

x = @SVector randn(4)
prob = ODEProblem(nelson,x,(0,1000.),2.0)

system = ContinuousDynamicalSystem(prob)
ds = CoupledODEs(system, diffeq)
##
u0s = [ SVector(p,1,0,0) for p ∈ LinRange(1,15,3) ]
# inputs
trs = [trajectory(ds, 10^5, u0)[1][:, SVector(1,2,3)] for u0 ∈ u0s]
j = 2 # the dimension of the plane

figure, ax3D, ax2D = interactive_poincaresos_scan(trs, j; linekw = (transparency = true,), scatterkw = (alpha=0.5,markersize=4,), colors=[:red,:green,:blue])

figure
##
trajectory(ds, 10000, u0s[1])[1][:,SVector(1,2,3)]