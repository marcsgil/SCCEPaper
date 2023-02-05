using SCCanonicalEnsemble
using FastGaussQuadrature

include("../plotting_functions.jl")
include("../ex_functions.jl")
##
H(x,χ) = χ*x[1]^2+(1-exp(-x[2]))^2/(4χ)

function f!(du,u,χ)
    du[1] = -4χ*u[3]
    du[2] = ( -exp(-u[4])*cos(0.5*u[1]) + exp(-2*u[4])*cos(u[1]) )/χ
    du[3] = (exp(-u[4])*sin(0.5*u[1]) - exp(-2*u[4])*sin(u[1]) )/(2χ)
    du[4] = -χ*u[2]
end

myexp(x) = exp(real(x))*cis(imag(x))
H′(x,χ) = χ*x[1]^2+(1-myexp(-x[2]))^2/(4χ)
g! = get_equations_of_motion(H′,1,"χ")

function getNodesAndWeights(χ,N=200)
    Ps,wPs = gausslegendre(N)
    Qs,wQs = gausschebyshev(N,3)

    [[√(1-Q^2)*P/(2χ),-log(1-Q)] for P in Ps, Q in Qs], [w1*w2 for w1 in wPs, w2 in wQs]
end
##
χ = .12
Es = [ (n+0.5)-χ*(n+0.5)^2 for n in 0:floor(Int64, 0.5*(1/χ-1))]
θ_min = .2
θ_max = 1
N = 32
θs_sc = LinRange(θ_min,θ_max,N)
θs_ex = LinRange(θ_min,θ_max,4N)
##
Us_ex = [ex_U(θ,Es) for θ in θs_ex]
Us_sc = solve_equations(θs_sc,χ,f!,getNodesAndWeights,H,output_func=energy_output,reduction=energy_reduction,callback=nothing)

@benchmark solve_equations(θ_max,χ,g!,getNodesAndWeights,H,callback=nothing)

comparison_plot(θs_ex,θs_sc,Us_ex,Us_sc)