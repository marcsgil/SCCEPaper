include("../../plot_config.jl")

θc(ϵ) = log(1 / √ϵ + √(1 / ϵ - 1)) / √(1 - ϵ)

ϵs = LinRange(0.001, 1, 1024)
f = Figure()
ax = Axis(f[1, 1], xticks=0:0.2:1, xlabel=L"\epsilon", ylabel=L"\omega s_c")
lines!(ax, ϵs, θc.(ϵs))
f
##
save("Plots/Morse/critical_time.pdf", f)
