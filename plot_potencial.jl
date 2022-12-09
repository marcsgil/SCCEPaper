using CairoMakie

V(q,χ) = χ*q^4/4+(1/2-χ)*q^2

qs = LinRange(-3,3,512)
##
f = Figure(fontsize=32)
ax = CairoMakie.Axis(f[1, 1],
    xlabel = L"q",
    xlabelsize=40,
    ylabel = L"v",
    ylabelsize=40
)
lines!(ax,qs,V.(qs,0),linewidth = 4,label=L"\chi=0")
lines!(ax,qs,V.(qs,.5),linewidth = 4,label=L"\chi=.5")
lines!(ax,qs,V.(qs,1),linewidth = 4,label=L"\chi=1")
lines!(ax,qs,V.(qs,2),linewidth = 4,label=L"\chi=2")
ylims!(ax,-1.5,5)
xlims!(ax,-3,3)
axislegend(position = :ct)
f