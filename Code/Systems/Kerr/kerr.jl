using SCCanonicalEnsemble
using ClassicalCanonicalEnsemble
using QuantumCanonicalEnsemble
using JLD2
##
H(x,χ) = sum(abs2,x)/2+χ*(sum(abs2,x)/2)^2
Es(χ) = [ (n+1/2) + χ*(n+1/2)^2 for n in 0:10^5 ]

sc_energy(θs::AbstractArray,χ) = energy_NF(θs,Polynomial([0,1,χ]),2000.)
sc_heat(θs::AbstractArray,χ) = heat_NF(θs,Polynomial([0,1,χ]),2000.)
##
θ_min = .2
θ_max = 10
N = 256
scatter_θs = LinRange(θ_min,θ_max,N)
line_θs = LinRange(θ_min,θ_max,N)
χ = .5
for χ ∈ (0.1,0.3,0.6,1.)
    quantum = quantum_energy(line_θs,χ,Es)
    sc = sc_energy(scatter_θs,χ) |> vec
    classical = classical_energy(line_θs,1,χ,hamiltonian=H)

    jldsave("Results/Kerr/energy_χ=$χ.jld2";scatter_θs,line_θs,χ,quantum,sc,classical)
end

χ = .5
for χ ∈ (0.1,0.3,0.6,1.)
    quantum = quantum_heat(line_θs,χ,Es)
    sc = sc_heat(scatter_θs,χ) |> vec
    classical = classical_heat(line_θs,1,χ,hamiltonian=H)

    jldsave("Results/Kerr/heat_χ=$χ.jld2";scatter_θs,line_θs,χ,quantum,sc,classical)
end