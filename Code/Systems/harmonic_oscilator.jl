using SCCanonicalEnsemble
##
H(x,ω) = x[1]^2/2 + x[2]^2/2

function f!(du,u,χ)
    du[1] = -2u[3]
    du[2] = -2u[4]
    du[3] = -u[1]/2
    du[4] = -u[2]/2
end
##
θ = .5

coth(θ/2)/2

energy_integrals(θ, 0, f!, H, [-10.,-10.], [10.,10.])
##
ub = ones(2) * 5
lb = - ub
energy_integrals2(θ, 0, f!, H, ub,lb)