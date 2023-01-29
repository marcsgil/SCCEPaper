using ProgressMeter

include("../sc_solution.jl")
include("../sc_functions.jl")
include("../ex_functions.jl")
include("../plotting_functions.jl")

V(q,μ) = (q[1]^2/2-q[2])^2 + μ*q[1]^2
H(x,μ) = (x[1]^2+x[2]^2)/2 + (x[3]^2/2-x[4])^2 + μ*x[3]^2
fy,fx = get_equations_of_motion(H,2)
##
μ = 1
xs = LinRange(-4.5,4.5,150)
ys = LinRange(-4,5,150)
Es,ψs = solveSchrodinger(xs,ys,V;nev=60,par=μ)

θ = 1
ex_U(θ,Es)

energyMonteCarlo(θ,μ,H,2,fy,fx,10^5)
##
θ_min = .5
θ_max = 3
N = 16
θs_sc = LinRange(θ_min,θ_max,N)
θs_ex = LinRange(θ_min,θ_max,4N)


Us_ex = [ex_U(θ,Es) for θ in θs_ex]
Us_sc = @showprogress [energyMonteCarlo(θ,μ,H,2,fy,fx,10^5) for θ in θs_sc]

comparison_plot(θs_ex,θs_sc,Us_ex,Us_sc,L"\mu = %$μ")