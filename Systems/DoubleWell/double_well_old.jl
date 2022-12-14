using FastGaussQuadrature, StaticArrays, ThreadsX, CairoMakie, Integrals
##
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
##
include("double_phase_space.jl")

V(q,χ) = χ*q^4/4+(1/2-χ)*q^2
H(x,χ) = x[1]^2/2 + V(x[2],χ)
fy(y,x,χ) = SA[-2x[1],-(χ*(4x[2]^2-3y[1]^2)/2+2(1-2χ))*x[2]]
fx(y,x,χ) = SA[(χ*(y[1]^2/4-3x[2]^2)-(1-2χ))*y[1]/2,-y[2]/2]

function quadrature_generator(θ,χ,N=20)
    coordinate_transformation((P,Q),θ,χ) = SA[√(2P/θ),(4Q/(θ*χ))^1/4]
    Ps,wPs = gausslaguerre(N,-1/2)
    Qs,wQs = gausslaguerre(N,-3/4)

    coordinate_transformation.(Iterators.product(Ps,Qs),θ,χ),
    map(prod,Iterators.product(wPs,wQs))
end

function Z_integrand(u,θ,χ,nodes,i)
    #=ifelse(exp(-θ*H(nodes[i],χ)) < 1e-3,0,
    exp(u.Δ- θ*(1/2-χ)*nodes[i][2]^2)*√abs(det(u.jac_x)))=#
    exp(u.Δ- θ*(1/2-χ)*nodes[i][2]^2)*√abs(det(u.jac_x))
end

function energy_output(sol,i,(nodes,weights),θ,par)
    replace([weights[i]*Z_integrand(sol[end],θ,χ,nodes,i),H(sol[end].x,par)],Inf=>0),false
end

function energy_reduction(sols,θ)
    sum(prod,eachcol(sols))/sum(first,view(sols,1,:))
    #sols
end

##
χ = 0.1
θ_min = .2
θ_max = 1
N = 16
θs_sc = LinRange(θ_min,θ_max,N)
θs_ex = LinRange(θ_min,θ_max,4N)
##
U_ex = EX_U(θs_ex,V,χ)
#U_sc = ThreadsX.map(θ->calculate_expectation(θ,fy,fx,quadrature_generator,energy_output,energy_reduction,χ), θs_sc)
calculate_expectation(1,fy,fx,quadrature_generator,energy_output,energy_reduction,χ)
findall(isinf,Array(U_sc[1]))

#sum(prod,eachcol(U_sc[1]))/sum(first,view(U_sc[1],1,:))
sum(first,view(U_sc[1],2,:))
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
axislegend()
px = .45
py = .91
text!(ax,px*θs_ex[end]+(1-px)θs_ex[1],py*U_sc[1]+(1-py)U_sc[end],text=L"\chi=%$χ",textsize=36)
f
##
function quadrature_generator(θ,χ,Ps,Qs)
    coordinate_transformation((P,Q)) = SA[P,Q]
    #Ps = LinRange(-5,5,128)
    #Qs = LinRange(-5,5,128)

    coordinate_transformation.(Iterators.product(Ps,Qs)),1
end

function Z_integrand(u,θ,χ,nodes,i)
    prefactor = exp(-θ*H(nodes[i],χ))
    ifelse(prefactor < 1e-3,0,prefactor*exp(u.Δ)*√abs(det(u.jac_x)))
end

function plot_output(sol,i,(nodes,weights),θ,par)
    Z_integrand(sol[end],θ,χ,nodes,i),false
end

function plot_reduction(sols,θ)
    reshape(sols,128,128)
end
##
χ = 0.5
Ps = LinRange(-200,200,128)
Qs = LinRange(-200,200,128)

Zs = calculate_expectation(.2,fy,fx,(θ,χ)->quadrature_generator(θ,χ,Ps,Qs),plot_output,plot_reduction,χ)

heatmap(Zs,colormap=:hot,axis=(;aspect=1),figure=(;resolution=(500,500)))

Zs

nodes,_ = quadrature_generator(.2,χ)

nodes

findall(iszero,Zs)