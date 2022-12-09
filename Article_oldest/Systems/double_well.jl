using CairoMakie
using Integrals,IntegralsCubature,IntegralsCuba

include("../thermo_quantities.jl")
include("../double_phase_space.jl")
include("../solveSchrodinger.jl")

V(q,depth) = depth*(q^2-1)^2
H(x,depth) = x[1]^2/depth + V(x[2],depth)
fy(y,x,depth) = [ -4x[1]/depth, 2depth*x[2]*(3*y[1]^2+4*(1-x[2])^2) ]
fx(y,x,depth) = [ 2depth*y[1]*(y[1]^2/4-3x[2]^2+1),-y[2]/depth ]

function EX_U(θ::Number,depth)
    xs = LinRange(-10,10,2048)
    Es,ψs = solveSchrodinger(xs,V;mass=depth/2,par=depth)
    sum(E->E*exp(-θ*E),Es)/sum(E->exp(-θ*E),Es)
end

function EX_U(θs::AbstractArray,depth)
    xs = LinRange(-10,10,2048)
    Es,ψs = solveSchrodinger(xs,V;mass=depth/2,par=depth)
    map(θ->sum(E->E*exp(-θ*E),Es)/sum(E->exp(-θ*E),Es),θs)
end

SC_U(θ::Number,depth) = energy(θ::Number,depth,fy,fx,H)
SC_U(θs::AbstractArray,depth) = energy(θs,depth,fy,fx,H)

function CL_energy_integrand!(dX,X,(θ,depth))
    Threads.@threads for n in axes(X,2)
        dX[:,n] = exp(-θ*H(view(X,:,n),depth))*[1,H(view(X,:,n),depth)]
    end
end

function CL_U(θ,depth)
    prob = IntegralProblem(CL_energy_integrand!,-20*ones(2),20*ones(2),(θ,depth),nout=2)
    sol = solve(prob,CubatureJLh(),reltol=1e-3,abstol=1e-4)
    sol[2]/sol[1]
end
##
depth=.01
EX_U(.5,depth)
SC_U(.5,depth)
CL_U(.5,depth)
θs = LinRange(.02,.3,50)

EX_Us = EX_U(θs,depth)
SC_Us = SC_U(θs,depth)
Cl_Us = CL_U.(θs,depth)
##
fig,ax,l = lines(θs,EX_Us,label="Exact")
lines!(ax,θs,SC_Us,label="SC")
lines!(ax,θs,Cl_Us,label="Classical")
ylims!(ax,1,20.05)
axislegend()
fig
##
error = map((ap,ex)->abs(1-ap/ex),EX_Us,SC_Us)
fig,ax,l = lines(θs,error)
#ylims!(ax,0,.08)
fig
