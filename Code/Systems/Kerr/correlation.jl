include("../../plot_config.jl")
E(n, χ) = (n + 1 / 2) + χ * (n + 1 / 2)^2

boltzmann_factor(n,β,χ) = exp(-β * E(n, χ) )


function position_correlation(t, β, χ, N)
    (sum(n -> exp(-β * E(n, χ)) * ((n + 1) * cis(t * (E(n + 1, χ) - E(n, χ))) + n * cis(t * (E(n - 1, χ) - E(n, χ)))), 0:N )
    / (2 * sum(n -> exp(-β * E(n, χ)), 0:N)))
end

function position_correlation(t, β, χ, N)
    coth(β/2)*cos(t)/2 + im *sin(t)/2
end
##
ts = LinRange(0,2π,128)
cs = [position_correlation(t, 1, .0, 10^4) for t ∈ ts]

fig,ax = lines(ts,real.(cs))
lines!(ax,ts,imag.(cs))
fig