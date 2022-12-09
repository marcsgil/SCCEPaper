using SemiClassicalCanonicalEnsemble
using CairoMakie
using StaticArrays

include("solveSchrodinger.jl")

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

V(q,χ) = χ*q^4/4+(1/2-χ)*q^2
H(x,χ) = x[1]^2/2 + V(x[2],χ)
fy(y,x,χ) = SA[-2x[1],-(χ*(4x[2]^2-3y[1]^2)/2+2(1-2χ))*x[2]]
fx(y,x,χ) = SA[(χ*(y[1]^2/4-3x[2]^2)-(1-2χ))*y[1]/2,-y[2]/2]
fy(y,x,ω) = -2*x
fx(y,x,ω) = -y/2
##
χ = 0.3
θ_min = .5
θ_max = 1
N = 16
θs_sc = LinRange(θ_min,θ_max,N)
θs_ex = LinRange(θ_min,θ_max,4N)
##
#U_sc = energy(θs_sc[1],χ,fy,fx,H) 
U_sc = energy(θs_sc,χ,fy,fx,H)
U_ex = EX_U(θs_ex,V,χ)
##
f = Figure(fontsize=24)
ax = Axis(f[1, 1],
    xlabel = L"θ",
    xlabelsize=32,
    ylabel = "Energy",
    ylabelsize=24
)
lines!(ax,θs_ex,U_ex,label="Exact",color=:black)
scatter!(ax,θs_sc,U_sc,label="SC",color=:red,marker=:diamond)
axislegend()
px = .45
py = .91
text!(ax,px*θs_ex[end]+(1-px)θs_ex[1],py*U_sc[1]+(1-py)U_sc[end],text=L"\chi=%$χ",textsize=36)
f