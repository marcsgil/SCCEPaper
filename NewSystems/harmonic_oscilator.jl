using StaticArrays

include("../sc_solution.jl")
include("../sc_functions.jl")

H(x,ω) = ω*sum(x->x^2,x)/2
fy,fx =  get_equations_of_motion(H,1)

exact_U(θ,ω) = coth(θ*ω/2)*ω/2
θ,ω = rand(),rand()
energyMonteCarlo(θ,ω,H,1,fy,fx,10^5)

exact_U(θ,ω)