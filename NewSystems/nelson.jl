include("../sc_solution.jl")
include("../sc_functions.jl")
include("../ex_functions.jl")
include("../plotting_functions.jl")

V(q,μ) = (q[1]^2/2-q[2])^2 + μ*q[1]^2
H(x,μ) = (x[1]^2+x[2]^2)/2 + (x[3]^2/2-x[4])^2 + μ*x[3]^2
