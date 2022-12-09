using Integrals
using CairoMakie

f(t,R) = R*cis(.5*R^2*cis(2t))

function I(R)
    prob = IntegralProblem(f,0,π/4,R)
    solve(prob,HCubatureJL()).u
end

I(4*10^2) |> abs2
Rs = LinRange(0,100,100)

modules = abs2.(I.(Rs))

lines(Rs,modules)
