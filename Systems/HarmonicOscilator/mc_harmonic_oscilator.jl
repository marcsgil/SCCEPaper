using FastGaussQuadrature, StaticArrays

include("../../double_phase_space.jl")
include("../../sampling.jl")

H(x,ω) = ω*x⋅x/2
H2(x,ω) = H(x,ω)^2-ω^2/4
fy(y,x,ω) = -2ω*x
fx(y,x,ω) = -ω*y/2

#=function Z_integrand(u,θ)
    exp(u.Δ)*√abs(det(u.jac_x))
end

function energy_output(sol,i,(nodes,weights),θ,par)
    [Z_integrand(sol[end],θ),H(sol[end].x,par)],false
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
    sum(prod,eachcol(sols))/sum(first,view(sols,1,:))
end

function quadrature_generator(θ,ω,N=10^5)   
    Ps = rand(Normal(0,1/√(θ*ω)),N)
    Qs = sample(q-> exp(-θ*ω*q^2/2),N)
    [SA[Ps[n],Qs[n]] for n in 1:N], ones(N)
end
##
exact_U(θ,ω) = coth(θ*ω/2)*ω/2
θ,ω = rand(),rand()
test = solve_equations(θ,ω,fy,fx,quadrature_generator,output_func=energy_output,reduction=energy_reduction)

exact_U(θ,ω)
##
ω = 1
θ_min = .5
θ_max = 4
N = 16
θs_sc = LinRange(θ_min,θ_max,N)
θs_ex = LinRange(θ_min,θ_max,4N)
##
U_ex = exact_U.(θs_ex,ω)
U_sc = map(θ->solve_equations(θ,ω,fy,fx,quadrature_generator,output_func=energy_output,reduction=energy_reduction), θs_sc)
##
using CairoMakie
f = Figure(fontsize=24)
ax = CairoMakie.Axis(f[1, 1],
    xlabel = L"θ",
    xlabelsize=32,
    ylabel = "Energy",
    ylabelsize=24
)
lines!(ax,θs_ex,U_ex,label="Exact",color=:black)
scatter!(ax,θs_sc,U_sc,label="SC",color=:red,marker=:diamond)
#lines!(ax,θs_sc,U_sc,label="SC",color=:red)
#lines!(ax,θs_ex,U_cl,label="Classical",color=:blue)
axislegend()
px = .45
py = .91
#text!(ax,px*θs_ex[end]+(1-px)θs_ex[1],py*U_ex[1]+(1-py)U_ex[end],text=L"\chi=%$χ",textsize=36)
#ylims!(ax,5U_ex[end],1.1*U_ex[1])
f