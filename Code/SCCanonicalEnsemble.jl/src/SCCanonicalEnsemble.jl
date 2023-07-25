module SCCanonicalEnsemble

using Reexport

@reexport using DifferentialEquations
using SparseDiffTools,LinearAlgebra,Integrals,IntegralsCubature,ProgressMeter,Tullio

#=include("sc_functions.jl")
export energy_output,energy_outputMonteCarlo,energy_reduction,analysis_output,
caustic_callback,energyMonteCarlo,heat_reduction,heat_output,heatMonteCarlo=#

include("double_hamiltonian_prob.jl")
export integrand

regularize(x) = isfinite(sum(x)) ? x : zero(x)

include("integrals.jl")
export energy_integrals, heat_integrals, energy_integrals2

include("projections.jl")
export wigner_output

include("quadrature.jl")
export energy_quadrature, heat_quadrature, caustic_callback, strong_callback

using OnlineStats, AdvancedMH, Distributions, MCMCChains, Statistics
include("monte_carlo.jl")
export energy_mc, heat_mc

using Symbolics
include("symbolic_eqs_of_motion.jl")
export get_equations_of_motion,squared_hamiltonian_symbol

@reexport using Polynomials, StaticArrays
include("normal_forms.jl")
export energy_NF,heat_NF


end
