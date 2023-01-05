include("../sc_solution.jl")
include("../sc_functions.jl")
include("../ex_functions.jl")
include("../plotting_functions.jl")

V(q,χ) = χ*q^4/4+(1/2-χ)*q^2
H(x,χ) = x[1]^2/2 + V(x[2],χ)
fy,fx = get_equations_of_motion(H,1)
##
χ = .5
xs = LinRange(-10,10,2048)
Es,ψs = solveSchrodinger(xs,V;par=χ)

θ = 2
ex_U(θ,Es)

energyMonteCarlo(θ,χ,H,1,fy,fx,10^5)
##
χ = 1
θ_min = .2
θ_max = 4
N = 16
θs_sc = LinRange(θ_min,θ_max,N)
θs_ex = LinRange(θ_min,θ_max,4N)

xs = LinRange(-10,10,2048)
Es,ψs = solveSchrodinger(xs,V;par=χ)

Us_ex = [ex_U(θ,Es) for θ in θs_ex]
Us_sc = [energyMonteCarlo(θ,χ,H,1,fy,fx,10^5) for θ in θs_sc]

comparison_plot(θs_ex,θs_sc,Us_ex,Us_sc)