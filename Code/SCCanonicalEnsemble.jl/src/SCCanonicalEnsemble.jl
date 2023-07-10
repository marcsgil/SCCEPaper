module SCCanonicalEnsemble

using Reexport

@reexport using DifferentialEquations
using SparseDiffTools,LinearAlgebra,Integrals,IntegralsCubature

#=include("sc_functions.jl")
export energy_output,energy_outputMonteCarlo,energy_reduction,analysis_output,
caustic_callback,energyMonteCarlo,heat_reduction,heat_output,heatMonteCarlo=#

include("double_hamiltonian_prob.jl")

include("integrals.jl")
export energy_integrals, heat_integrals

using Symbolics
include("symbolic_eqs_of_motion.jl")
export get_equations_of_motion,squared_hamiltonian_symbol

include("quadrature.jl")
export energy_quadrature

@reexport using Polynomials, StaticArrays
include("normal_forms.jl")
export energy_NF,heat_NF


end
