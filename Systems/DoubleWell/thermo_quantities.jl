using StaticArrays, CairoMakie, Integrals, IntegralsCubature, Distributions

include("../../solveSchrodinger.jl")
include("../../double_phase_space.jl")
include("../../metropolisHastings.jl")
include("double_well.jl")

function EX_U(θ::Number,V,χ)
    xs = LinRange(-10,10,2048)
    Es,ψs = solveSchrodinger(xs,V;par=χ)
    sum(E->E*exp(-θ*E),Es)/sum(E->exp(-θ*E),Es)
end

function EX_U(θs::AbstractArray,V,χ)
    xs = LinRange(-10,10,2048)
    Es,ψs = solveSchrodinger(xs,V;par=χ)
    map(θ->sum(E->E*exp(-θ*E),Es)/sum(E->exp(-θ*E),Es),θs)
end

function quadrature_generator(θ,χ,N=10^5)   
    Ps = rand(Normal(0,1/√θ),N)
    Qs = metropolisHastings(q-> exp(-θ*V(q,χ)),N)
    [SA[Ps[n],Qs[n]] for n in 1:N], ones(N)
end

function Z_integrand(u,θ,χ,nodes,i)
    exp(u.Δ)*√abs(det(u.jac_x))
end

function energy_output(sol,i,θ,par,node,weight,)
    [weight*Z_integrand(sol[end],θ,χ,node,i),H(sol[end].x,par)],false
end

function energy_reduction(sols,θ)
    sum(prod,eachcol(sols))/sum(first,view(sols,1,:))
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
θ=3
χ=.5
U_ex = EX_U(θ,V,χ)
U_sc = solve_equations(θ,χ,fy,fx,(θ,χ)->quadrature_generator(θ,χ,10^6) ,output_func=energy_output,reduction=energy_reduction)
CL_U(θ,χ)
##
χ = 1
θ_min = .2
θ_max = 3
N = 16
θs_sc = LinRange(θ_min,θ_max,N)
θs_ex = LinRange(θ_min,θ_max,4N)
##
U_ex = EX_U(θs_ex,V,χ)
U_sc = map(θ->solve_equations(θ,χ,fy,fx,quadrature_generator,output_func=energy_output,reduction=energy_reduction), θs_sc)
U_cl = map(θ-> CL_U(θ,χ),θs_ex)
##
CairoMakie.activate!()
f = CairoMakie.Figure(fontsize=24)
ax = CairoMakie.Axis(f[1, 1],
    xlabel = L"θ",
    xlabelsize=32,
    ylabel = "Energy",
    ylabelsize=24
)
lines!(ax,θs_ex,U_ex,label="Exact",color=:black)
scatter!(ax,θs_sc,U_sc,label="SC",color=:red,marker=:diamond)
#lines!(ax,θs_sc,U_sc,label="SC",color=:red)
lines!(ax,θs_ex,U_cl,label="Classical",color=:blue)
axislegend()
px = .45
py = .91
text!(ax,px*θs_ex[end]+(1-px)θs_ex[1],py*U_ex[1]+(1-py)U_ex[end],text=L"\chi=%$χ",textsize=36)
#ylims!(ax,5U_ex[end],1.1*U_ex[1])
f