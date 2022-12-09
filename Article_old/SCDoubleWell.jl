include("double_phase_space.jl")
using FastGaussQuadrature,StaticArrays

V(q,χ) = χ*q^4/4+(1/2-χ)*q^2
H(x,χ) = x[1]^2/2 + V(x[2],χ)
fy(y,x,χ) = SA[-2x[1],-(χ*(4x[2]^2-3y[1]^2)/2+2(1-2χ))*x[2]]
fx(y,x,χ) = SA[(χ*(y[1]^2/4-3x[2]^2)-(1-2χ))*y[1]/2,-y[2]/2]

coord_transformation((P,Q),θ,χ) = SA[√(2P/θ),(4Q/θ*χ)^1/4]

function get_nws(N,θ,χ)
    Ps,wPs = gausslaguerre(N,-1/2)
    Qs,wQs = gausslaguerre(N,-3/4)

    coord_transformation.(Iterators.product(Ps,Qs),θ,χ), [w1*w2 for w1 in wPs, w2 in wQs]
end

function energy_output(sol,i,θ,par,(nodes,weights),H)
    Z = weights[i]*exp( sol[end].Δ - θ*(1/2-χ)*nodes[i][2]^2 )*√abs(det(sol[end].jac_x))
    [Z,Z*H(sol[end].x,par)],false
end

function energy_reduction(sol)
    sum( n-> isinf(sol[1,n]*sol[2,n]) ? 0 : sol[1,n]*sol[2,n],1:length(sol))/sum( n->isinf(sol[1,n]) ? 0 : sol[1,n],1:length(sol))
end

θ = .1
χ = .5
nodes,weights = get_nws(100,θ,χ)
test = calculate_expectation(θ,fy,fx,(nodes,weights),energy_output,energy_reduction,χ)

test[1,2]

test |> Array
sum( x->x[1],test)
energy_reduction(test)

function U(θ,χ,N)
    ns,ws = get_nws(N,θ,χ)
end