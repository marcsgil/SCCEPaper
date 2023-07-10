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

energyMonteCarlo(θ,0,H,1,f!,10^6,callback=nothing)
##