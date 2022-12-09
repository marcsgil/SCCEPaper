module SemiClassicalCanonicalEnsemble


include("equations_of_motion.jl")


include("double_phase_space.jl")

using Integrals,IntegralsCubature,IntegralsCuba,ThreadsX
include("thermo_quantities.jl")

test() = "Success!"
test2() = "Another Success!"

export energy

end
