module SCCanonicalEnsemble

using Reexport

@reexport using DifferentialEquations
using SparseDiffTools,LinearAlgebra,MetropolisHastings

include("sc_solution.jl")
export caustic_callback,strong_callback,disc_caustic_callback,
solve_equations,extract_jac_x,phase_space_dim,extract_det_jac,annul!

include("sc_functions.jl")
export energy_output,energy_outputMonteCarlo,energy_reduction,analysis_output,
caustic_callback,energyMonteCarlo,heat_reduction,heat_output,heatMonteCarlo

using Symbolics
include("symbolic_eqs_of_motion.jl")
export get_equations_of_motion,squared_hamiltonian_symbol

using StaticArrays
@reexport using Polynomials
include("normal_forms.jl")
export energy_NF,heat_NF


end
