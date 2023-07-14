include("../../plot_config.jl")

V(q) = (1-exp(-q))^2

qs = LinRange(-.8,4,1024)
f = Figure()
ax = Axis(f[1, 1], xticks=0:0.2:1, xlabel=L"q", ylabel=L"V(q)", 
yticks = (0:0.5:1.5, [L"0", L"D/2", L"D", L"3D/2"]))
lines!(ax, ϵs, V.(qs))
f
##
save("Plots/Morse/morse_potential.pdf", f)