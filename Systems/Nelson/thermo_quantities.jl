using StaticArrays, CairoMakie, Integrals, IntegralsCubature, Distributions

include("../../solveSchrodinger.jl")
include("../../double_phase_space.jl")
include("../../metropolisHastings.jl")
include("nelson.jl")

function quadrature_generator(θ,μ,N=10^5)
    xs = metropolisHastings(x -> exp(-θ*H(x,μ) ),N,4)
    [SVector{4}(xs[:,n]) for n in 1:N], ones(N)
end

function Z_integrand(u,θ,μ,nodes,i)
    exp(u.Δ)*√abs(det(u.jac_x))
end

function energy_output(sol,i,θ,par,node,weight)
    [weight*Z_integrand(sol[end],θ,μ,node,i),H(sol[end].x,par)],false
end

function energy_reduction(sols,θ)
    sum(prod,eachcol(sols))/sum(first,view(sols,1,:))
end

function CL_energy_integrand!(dX,X,(θ,μ))
    Threads.@threads for n in axes(X,2)
        dX[:,n] = exp(-θ*H(view(X,:,n),μ))*[1,H(view(X,:,n),μ)]
    end
end

function CL_U(θ,μ)
    prob = IntegralProblem(CL_energy_integrand!,-20*ones(4),20*ones(4),(θ,μ),nout=2)
    sol = solve(prob,CubatureJLh(),reltol=1e-4,abstol=1e-5)
    sol[2]/sol[1]
end
##
θ=.1
μ= 2

U_sc = solve_equations(θ,μ,fy,fx,quadrature_generator,output_func=energy_output,reduction=energy_reduction)
CL_U(θ,μ)
##
nodes,weights = quadrature_generator(θ,μ,10^6)

sum(x->H(x,μ)*exp(-θ*H(x,μ)),nodes)/10^6
##
xs = LinRange(-4.5,4.5,128)
ys = LinRange(-3,4,128)
Es,ψs = solveSchrodinger(xs,ys,V;nev=50,par=μ)

Es
ψs
##

heatmap(xs,ys,abs2.(reshape(ψs[:,9],length(xs),length(ys))))
##
U(θ,Es) = sum(E->E*exp(-θ*E),Es)/sum(E->exp(-θ*E),Es)

U(θ,Es)