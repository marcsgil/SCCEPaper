module SCCanonicalEnsemble

using Reexport

@reexport using DifferentialEquations
using SparseDiffTools,LinearAlgebra

include("sc_solution.jl")
export caustic_callback,strong_callback,disc_caustic_callback,
solve_equations,extract_jac_x,phase_space_dim,extract_det_jac,annul!

include("metropolis_hastings.jl")

include("sc_functions.jl")
export energy_output,energy_outputMonteCarlo,energy_reduction,analysis_output,
caustic_callback,energyMonteCarlo

using Symbolics
include("symbolic_eqs_of_motion.jl")
export get_equations_of_motion

using StaticArrays
@reexport using Polynomials
include("normal_forms.jl")
export energy_NF,heat_NF


end
