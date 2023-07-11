using SCCanonicalEnsemble

V(q,χ) = χ*q^4/4+(1/2-χ)*q^2
H(x,χ) = x[1]^2/2 + V(x[2],χ)
f! = get_equations_of_motion(H,1,"χ")

energy_mc(1, 1, f!, H, 1,maxiters=2*10^6)