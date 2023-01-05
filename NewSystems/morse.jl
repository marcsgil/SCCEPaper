using FastGaussQuadrature

include("../sc_solution.jl")
include("../sc_functions.jl")
include("../ex_functions.jl")
include("../plotting_functions.jl")
##
myexp(x) = exp(real(x))*cis(imag(x))
H(x,χ) = χ*x[1]^2+(1-myexp(-x[2]))^2/(4χ)
fy,fx = get_equations_of_motion(H,1)

function Z_integrand(u,energy_term)
    exp( u.Δ - energy_term )*√abs(det(u.jac_x))
end

function energy_output(sol,i,θs,par,nodes,weights,H)
    output = Array{eltype(θs)}(undef,2,length(θs))

    E = H(nodes[i],par)
    for (n,u) in enumerate(sol.u)
        output[1,n] = weights[i]*Z_integrand(u,θs[n]*E)
        output[2,n] = H(u.x,par)
    end
    output,false
end

function energy_reduction(sols,θs)
    f(collection) = sum(prod,eachcol(collection))/sum(first,view(collection,1,:))
    map(n->f(view(sols,:,n,:)),axes(sols,2))
end

function getNodesAndWeights(χ,N)
    Ps,wPs = gausslegendre(N)
    Qs,wQs = gausschebyshev(N,3)

    [SA[√(1-Q^2)*P/(2χ),-log(1-Q)] for P in Ps, Q in Qs], [w1*w2 for w1 in wPs, w2 in wQs]
end
##
χ = .12
θ_min = .5
θ_max = 5
N = 16
θs_sc = LinRange(θ_min,θ_max,N)
θs_ex = LinRange(θ_min,θ_max,4N)

Es = [ (n+0.5)-χ*(n+0.5)^2 for n in 0:floor(Int64, 0.5*(1/χ-1))]

Us_ex = [ex_U(θ,Es) for θ in θs_ex]
Us_sc = solve_equations(θs_sc,χ,fy,fx,χ -> getNodesAndWeights(χ,200),H,
output_func=energy_output,reduction=energy_reduction,callback=weak_callback)

comparison_plot(θs_ex,θs_sc,Us_ex,Us_sc)