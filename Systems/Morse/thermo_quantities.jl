using StaticArrays, FastGaussQuadrature, CairoMakie

include("../../double_phase_space.jl")
include("morse.jl")

function EX_U(θ::Number,χ)
    E(n,χ) = (n+0.5)-χ*(n+0.5)^2
    Nmax = floor(Int64, 0.5*(1/χ-1))
    sum(n->E(n,χ)*exp(-θ*E(n,χ)),0:Nmax)/sum(n->exp(-θ*E(n,χ)),0:Nmax)
end

EX_U(θs::AbstractArray,χ) = EX_U.(θs,χ)

function getNodesAndWeights(χ,N=100)
    Ps,wPs = gausslegendre(N)
    Qs,wQs = gausschebyshev(N,3)

    [SA[√(1-Q^2)*P/(2χ),-log(1-Q)] for P in Ps, Q in Qs], [w1*w2 for w1 in wPs, w2 in wQs]
end

function Z_integrand(u,energy_term)
    exp( u.Δ - energy_term )*√abs(det(u.jac_x))
end

function energy_output(sol,i,θs,par,node,weight)
    output = fill( ntuple(i->zero(eltype(θs)),2), length(θs) )
    E = H(node,par)
    for (n,u) in enumerate(sol.u)
        output[n] = (weight*Z_integrand(u,θs[n]*E),H(u.x,par))
    end
    output,false
end

function energy_reduction(sols,θs)
    f(collection) = sum( x->x[1]*x[2],collection)/sum( x->x[1],collection)
    map(n->f(view(sols,n,:)),eachindex(θs))
end

function CL_energy_integrand!(dX,X,(θ,χ))
    Threads.@threads for n in axes(X,2)
        dX[:,n] = exp(-θ*H(view(X,:,n),χ))*[1,H(view(X,:,n),χ)]
    end
end

function CL_U(θ,χ)
    prob = IntegralProblem(CL_energy_integrand!,-20*ones(2),20*ones(2),(θ,χ),nout=2)
    sol = solve(prob,CubatureJLh(),reltol=1e-4,abstol=1e-5)
    sol[2]/sol[1]
end
##
nplot = 32
θmin = .1
θmax = 6
θs_sc=LinRange(θmin,θmax,nplot)
θs_ex=LinRange(θmin,θmax,4nplot)
χ=.12
U_ex = EX_U(θs_ex,χ)
U_sc = solve_equations(θs_sc,χ,fy,fx,χ->getNodesAndWeights(χ,200),output_func=energy_output,reduction=energy_reduction)
#U_cl = map(θ-> CL_U(θ,χ),θs_ex)
##
f = Figure(fontsize=24)
ax = CairoMakie.Axis(f[1, 1],
    xlabel = L"θ",
    xlabelsize=32,
    ylabel = "Energy",
    ylabelsize=24
)
lines!(ax,θs_ex,U_ex,label="Exact",color=:black)
scatter!(ax,θs_sc,U_sc,label="SC",color=:red,marker=:diamond)
#lines!(ax,θs,U_sc,label="SC",color=:red)
#lines!(ax,θs_ex,U_cl,label="Classical",color=:blue)
axislegend()
px = .45
py = .91
text!(ax,px*θs_ex[end]+(1-px)θs_ex[1],py*U_ex[1]+(1-py)U_ex[end],text=L"\chi=%$χ",textsize=36)
#ylims!(ax,5U_ex[end],1.1*U_ex[1])
f