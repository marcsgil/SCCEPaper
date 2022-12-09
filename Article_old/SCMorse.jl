module SCMorse

using StaticArrays,FastGaussQuadrature,ThreadedIterables

include("../../CompleteSemiclassicalApprox/double_phase_space.jl")
include("../../CompleteSemiclassicalApprox/thermo_quantities.jl")
include("../../utils.jl")

fy(y,x,χ) = SA[-4χ*x[1], ( -exp(-x[2])*cos(0.5*y[1]) + exp(-2*x[2])*cos(y[1]) )/χ]
fx(y,x,χ) = SA[ (exp(-x[2])*sin(0.5*y[1]) - exp(-2*x[2])*sin(y[1]) )/(2χ), -χ*y[2]]
H(x,χ) = χ*x[1]^2+(1-exp(-x[2]))^2/(4χ)
H2(x,χ) = H(x,χ)^2 - (2exp(-2x[2])-exp(-x[2]))/4
coord_transformation((P,Q),χ) = SA[√(1-Q^2)*P/(2χ),-log(1-Q)]

function U(θs::AbstractArray,χ::AbstractFloat,(nodes,weights))
    coord_transformation((P,Q),χ) = SA[√(1-Q^2)*P/(2χ),-log(1-Q)]
    transformed_nodes = map(node->coord_transformation(node,χ),nodes)

    morse_energy_output(sol,i,(transformed_nodes,weights),θs,par) = energy_output(sol,i,(transformed_nodes,weights),θs,par,H)

    calculate_expectation(θs::AbstractArray,fy,fx,(transformed_nodes,weights),morse_energy_output,energy_reduction,χ)
end

U(θs::AbstractArray,χ::AbstractFloat,N::Int) = U(θs,χ,twoD_nws(get_last_half.(gausslegendre(N)), gausschebyshev(N,3)))

function U(θs::AbstractArray,χs::AbstractArray,N::Int)
    nws = twoD_nws(get_last_half.(gausslegendre(N)), gausschebyshev(N,3))
    result = tmap(χ->U(θs,χ,nws),χs)
    reduce(hcat,result)
end

function C(θs::AbstractArray,χ::AbstractFloat,(nodes,weights))
    coord_transformation((P,Q),χ) = SA[√(1-Q^2)*P/(2χ),-log(1-Q)]
    transformed_nodes = map(node->coord_transformation(node,χ),nodes)

    morse_heat_output(sol,i,(transformed_nodes,weights),θs,par) = heat_output(sol,i,(transformed_nodes,weights),θs,par,H,H2)

    calculate_expectation(θs::AbstractArray,fy,fx,(transformed_nodes,weights),morse_heat_output,heat_reduction,χ)
end

C(θs::AbstractArray,χ::AbstractFloat,N::Int) = C(θs,χ,twoD_nws(get_last_half.(gausslegendre(N)), gausschebyshev(N,3)))

function C(θs::AbstractArray,χs::AbstractArray,N::Int)
    nws = twoD_nws(get_last_half.(gausslegendre(N)), gausschebyshev(N,3))
    result = tmap(χ->C(θs,χ,nws),χs)
    reduce(hcat,result)
end

U([1e-4],1e-4,3)
C([1e-4],1e-4,3)

end