using SCCanonicalEnsemble
using FastGaussQuadrature,StaticArrays

include("../plotting_functions.jl")
include("../ex_functions.jl")
##
myexp(x) = exp(real(x))*cis(imag(x))
H(x,χ) = χ*x[1]^2+(1-myexp(-x[2]))^2/(4χ)
#H(x,χ) = χ*x[1]^2+(1-exp(-x[2]))^2/(4χ)
fy,fx = get_equations_of_motion(H,1)
gy(y,x,χ) = SA[-4χ*x[1], ( -exp(-x[2])*cos(0.5*y[1]) + exp(-2*x[2])*cos(y[1]) )/χ]
gx(y,x,χ) = SA[ (exp(-x[2])*sin(0.5*y[1]) - exp(-2*x[2])*sin(y[1]) )/(2χ), -χ*y[2]]

function g(u,χ,t)
    SA[-4χ*u[3],
    ( -exp(-u[4])*cos(0.5*u[1]) + exp(-2*u[4])*cos(u[1]) )/χ,
    (exp(-u[4])*sin(0.5*u[1]) - exp(-2*u[4])*sin(u[1]) )/(2χ),
    -χ*u[2],
    (exp(-u[4])*sin(0.5*u[1]) - exp(-2*u[4])*sin(u[1]) )/(2χ)*u[1] - χ*u[2]*u[2]]
end

nodes,weights = getNodesAndWeights(χ,200)

X = nodes[1]
u0 = SCCanonicalEnsemble.u02(nodes[1])


f = get_equations_of_motion2(H,1)

function getNodesAndWeights(χ,N)
    Ps,wPs = gausslegendre(N)
    Qs,wQs = gausschebyshev(N,3)

    [SA[√(1-Q^2)*P/(2χ),-log(1-Q)] for P in Ps, Q in Qs], [w1*w2 for w1 in wPs, w2 in wQs]
end
##
χ = .12
θ_min = .1
θ_max = 5
N = 32
θs_sc = LinRange(θ_min,θ_max,N)
θs_ex = LinRange(θ_min,θ_max,4N)
##
@benchmark solve_equations(θ_max,χ,gy,gx,χ -> getNodesAndWeights(χ,200),H,callback=nothing)
@benchmark solve_equations2(θ_max,χ,g,χ -> getNodesAndWeights(χ,200),H,callback=nothing)
##
Es = [ (n+0.5)-χ*(n+0.5)^2 for n in 0:floor(Int64, 0.5*(1/χ-1))]

Us_ex = [ex_U(θ,Es) for θ in θs_ex]
Us_sc = solve_equations(θs_sc,χ,fy,fx,χ -> getNodesAndWeights(χ,200),H,
output_func=energy_output,reduction=energy_reduction,callback=caustic_callback)

comparison_plot(θs_ex,θs_sc,Us_ex,Us_sc)
##

getNodesAndWeights(N) = [ SA[-.01,-log(2)+.01],SA[.01,-log(2)+.01] ],nothing
##
sols = solve_equations(2.04,χ,fy,fx,getNodesAndWeights,H,callback=nothing)

using LinearAlgebra


eigvals(sols[1].u[1].jac_x)
eigvals(sols[2].u[1].jac_x )
det(sols[1].u[1].jac_x )
det(sols[2].u[1].jac_x )

eigvals(sols[1].u[1].jac_y*inv(sols[1].u[1].jac_x))
eigvals(sols[2].u[1].jac_y*inv(sols[2].u[1].jac_x))
sols[1].u[1].jac_y*inv(sols[1].u[1].jac_x)
sols[2].u[1].jac_y*inv(sols[2].u[1].jac_x)

sqrt(complex(sols[1].u[1].jac_y |> det))/√(abs(sols[1].u[1].jac_y |> det))
sqrt(complex(sols[2].u[1].jac_y |> det))/√(abs(sols[2].u[1].jac_y |> det))