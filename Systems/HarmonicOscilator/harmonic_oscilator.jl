using FastGaussQuadrature, StaticArrays

include("double_phase_space.jl")

H(x,ω) = ω*x⋅x/2
H2(x,ω) = H(x,ω)^2-ω^2/4
fy(y,x,ω) = -2ω*x
fx(y,x,ω) = -ω*y/2

function Z_integrand(u,θ)
    exp(u.Δ)*√abs(det(u.jac_x))
end

function quadrature_generator(θ,ω,N=20)
    coordinate_transformation((P,Q)) = SA[√(2P/(ω*θ)),√(2Q/(ω*θ))]
    nodes,weights = gausslaguerre(N,-1/2)
    coordinate_transformation.(Iterators.product(nodes,nodes)),
    map(w->w[1]*w[2],Iterators.product(weights,weights))
end

function energy_output(sol,i,(nodes,weights),θ,par)
    [weights[i]*Z_integrand(sol[end],θ),H(sol[end].x,par)],false
end

function energy_reduction(sols,θ)
    sum(prod,eachcol(sols))/sum(first,view(sols,1,:))
end
##
θ,ω = rand(),rand()
test = calculate_expectation(θ,fy,fx,quadrature_generator,energy_output,energy_reduction,ω)

exact_U(θ,ω) = coth(θ*ω/2)*ω/2
exact_U(θ,ω)