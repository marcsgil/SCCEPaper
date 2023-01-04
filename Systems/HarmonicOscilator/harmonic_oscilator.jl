using FastGaussQuadrature, StaticArrays

include("../../double_phase_space.jl")

H(x,ω) = ω*x⋅x/2
H2(x,ω) = H(x,ω)^2-ω^2/4
fy(y,x,ω) = -2ω*x
fx(y,x,ω) = -ω*y/2

function quadrature_generator(θ,ω,N=20)
    coordinate_transformation((P,Q)) = SA[√(2P/(ω*θ)),√(2Q/(ω*θ))]
    nodes,weights = gausslaguerre(N,-1/2)
    coordinate_transformation.(Iterators.product(nodes,nodes)),
    map(w->w[1]*w[2],Iterators.product(weights,weights))
end

#=function Z_integrand(u,θ)
    exp(u.Δ)*√abs(det(u.jac_x))
end

function energy_output(sol,i,(nodes,weights),θ,par)
    [weights[i]*Z_integrand(sol[end],θ),H(sol[end].x,par)],false
end

function energy_reduction(sols,θ)
    sum(prod,eachcol(sols))/sum(first,view(sols,1,:))
end=#

function Z_integrand(u,θ,par,nodes,i)
    exp(u.Δ)*√abs(det(u.jac_x))
end

function energy_output(sol,i,node,weight,θ,par)
    [weight*Z_integrand(sol[end],θ,par,node,i),H(sol[end].x,par)],false
end

function energy_reduction(sols,θ)
    sum(prod,eachcol(sols))/sum(view(sols,1,:))
end
##
θ,ω = 1.5,1
test = solve_equations(θ,ω,fy,fx,quadrature_generator,output_func=energy_output,reduction=energy_reduction)

sum(prod,eachcol(test))/sum(view(test,1,:))

exact_U(θ,ω) = coth(θ*ω/2)*ω/2
exact_U(θ,ω)
##