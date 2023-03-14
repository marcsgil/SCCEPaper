using SCCanonicalEnsemble

H(x,χ) = χ*x[1]^2+(1-exp(-x[2]))^2/(4χ)
H2(x,χ) = H(x,χ)^2 - (2exp(-2x[2])-exp(-x[2]))/4

H3 = squared_hamiltonian_symbol(H,1)

H2([.1,.1],10)
H3([.1,.1],10)
