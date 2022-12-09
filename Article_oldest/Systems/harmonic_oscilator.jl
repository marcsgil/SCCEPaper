#include("../thermo_quantities.jl")
#include("../double_phase_space.jl")
include("../equations_of_motion.jl")
##
using LinearAlgebra
H(x,ω) = ω*x⋅x/2
symbolic_fy,symbolic_fx = equations_of_motion((p,q,ω) -> ω*(p^2+q^2)/2)

H2(x,ω) = H(x,ω)^2-ω^2/4
fy(y,x,ω) = -2ω*x
fx(y,x,ω) = -ω*y/2

fy([π,1/π],[exp(-1),exp(1)],log(2))
fx([π,1/π],[exp(-1),exp(1)],log(2))
symbolic_fy([π,1/π],[exp(-1),exp(1)],log(2))
symbolic_fx([π,1/π],[exp(-1),exp(1)],log(2))

using BenchmarkTools
@benchmark fy([π,1/π],[exp(-1),exp(1)],log(2))
@benchmark symbolic_fy([π,1/π],[exp(-1),exp(1)],log(2))

@benchmark fx([π,1/π],[exp(-1),exp(1)],log(2))
@benchmark symbolic_fx([π,1/π],[exp(-1),exp(1)],log(2))

EX_U(θ,ω) = ω*coth(θ*ω/2)/2
EX_C(θ,ω) = ((θ*ω/2)*csch(θ*ω/2))^2

SC_U(θ::Number,ω) = energy(θ,ω,fy,fx,H)
SC_U(θs::AbstractArray,ω) = energy(θs::AbstractArray,ω,fy,fx,H)
SC_C(θ::Number,ω) = heat(θ,ω,fy,fx,H,H2)
SC_C(θs::AbstractArray,ω) = heat(θs::AbstractArray,ω,fy,fx,H,H2)
##