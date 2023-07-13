using TaylorIntegration, SCCanonicalEnsemble
include("../../plot_config.jl")
V(q, μ) = (q[1]^2 / 2 - q[2])^2 + μ * q[1]^2
H(x, μ) = (x[1]^2 + x[2]^2) / 2 + (x[3]^2 / 2 - x[4])^2 + μ * x[3]^2
f! = get_equations_of_motion(H, 2)
g!(dest,src,par,t) = f!(dest,src,par)
x0 = [0, 0, 1.,1.] #the initial condition
t0 = 0.0     #the initial time
tmax = 100.0 #final time of integration
set_variables("z",numvars=4)
tv, xv, λv = lyap_taylorinteg(g!, x0, t0, tmax, 28, 1e-20, (1.); maxsteps=2000000);
##
fig,ax,ln = lines(tv, λv[:,1])
lines!(ax,tv, λv[:,2])
lines!(ax,tv, λv[:,3])
fig