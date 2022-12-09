using CairoMakie

include("Systems/harmonic_oscilator.jl")

ω = 1
θ = 1

EX_U(θ,ω)
SC_U(θ,ω)
##
θs = LinRange(.1,2,100)
Us = SC_U(θs,ω)
##
fig,ax,l = lines(θs,Us)
lines!(ax,θs,θ->EX_U(θ,ω))
ylims!(ax,-.1,10)
fig
##
EX_C(θ,ω)
SC_C(θ,ω)
##
θs = LinRange(.1,4,50)
Cs = SC_C(θs,ω)
##
fig,ax,l = lines(θs,Cs)
lines!(ax,θs,θ->EX_C(θ,ω))
ylims!(ax,-.1,1.05)
fig