using SCCanonicalEnsemble,QuantumCanonicalEnsemble
using SolveSchrodinger
using ProgressMeter

#include("../plotting_functions.jl")
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
χ = 1.
xs = LinRange(-10,10,2048)
Es,ψs = solveSchrodinger(xs,V;par=χ)

θ = 4.
quantum_heat(θ,Es)

heat_integrals(θ, χ, f!, H, [-10.,-10.], [10.,10.])
heat_mc(θ, χ, f!, H, 1,maxiters=10^7,batchsize=2^15)