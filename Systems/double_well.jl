using SCCanonicalEnsemble
using SolveSchrodinger
using ProgressMeter

include("../plotting_functions.jl")
##
V(q,χ) = χ*q^4/4+(1/2-χ)*q^2
H(x,χ) = x[1]^2/2 + V(x[2],χ)
f! = get_equations_of_motion(H,1,"χ")

#=function f!(du,u,χ)
    du[1] = -2*u[3]
    du[2] = -4u[4]*(1/2 - χ) - χ*(4u[4]*(u[4]^2 - u[1]^2/4) - 2*u[4]*u[1]^2)/2
    du[3] = χ*(-2*u[1]*u[4]^2 - u[1]*(u[4]^2 - u[1]^2/4))/2 - u[1]*(1/2 - χ)
    du[4] = -u[2]/2
    nothing
end=#
##
χ = .5
xs = LinRange(-10,10,2048)
Es,ψs = solveSchrodinger(xs,V;par=χ)
ex_U(θ,Es) = sum(E->E*exp(-θ*E),Es)/sum(E->exp(-θ*E),Es)

θ = .5
ex_U(θ,Es)

energyMonteCarlo(θ,χ,H,1,f!,10^6,callback=strong_callback)
##
N = 128
ps = LinRange(-12,12,N)
qs = LinRange(-6,6,N)
xs = [[p,q] for p in ps, q in qs]

θs = LinRange(0,3,64)
χ= .5

custom_condition(u,t,integrator) = u[end]-1
custom_callback = ContinuousCallback(custom_condition,annul!,save_positions=(false,false))

sols = solve_equations(θs,χ,f!,(par)->(xs,nothing),H,
abstol=1e-12,reltol=1e-12,alg=Vern8(),output_func=analysis_output,callback=CallbackSet(custom_callback,disc_caustic_callback))
analysis_plot(sols,ps,qs,θs,χ)
ma
##
χ = .5
θ_min = .2
θ_max = 3
N = 16
θs_sc = LinRange(θ_min,θ_max,N)
θs_ex = LinRange(θ_min,θ_max,4N)

xs = LinRange(-10,10,2048)
Es,ψs = solveSchrodinger(xs,V;par=χ)

Us_ex = [ex_U(θ,Es) for θ in θs_ex]
Us_sc = @showprogress [energyMonteCarlo(θ,χ,H,1,f!,10^6,callback=CallbackSet(custom_callback,disc_caustic_callback)) for θ in θs_sc]

comparison_plot(θs_ex,θs_sc,Us_ex,Us_sc,L"\chi = %$χ")