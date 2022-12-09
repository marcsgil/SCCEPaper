using CairoMakie
##
include("fft_utils.jl")
include("solve1DSchrodinger.jl")
include("waveFunc2WigFunc.jl")
##
V(q) = q^2
qs = LinRange(-5,5,2048)
Es,ψs = solve(V,qs,1/2,1/2)
lines(qs,ψs[:, 21])
lines(qs,abs2.(ψs[:, 22]))
Es[1]
##
ps_grid = DFTGrid(2048,5,1)
ps = direct_grid(ps_grid)
##
Ws = calculate_many_wigners(ψs,qs,ps_grid,7)
##
round2(x) = round(x,digits=2)
E = [p^2/2 + V(q)/2 for q in qs, p in ps]
Ns = [1,3,5,7]
#lim = .8*maximum(Ws[:,:,n])
lim = .16
fig =Figure(fontsize=20,resolution=(1200,1100))
ax1 = Axis(fig[1,1],xlabel=L"q",ylabel=L"p",title=L"E_{%$(Ns[1]-1)}=%$(round2(Es[Ns[1]]))"
    ,xlabelsize=25,ylabelsize=25,titlesize=25,aspect=1)
ax2 = Axis(fig[1,2],xlabel=L"q",ylabel=L"p",title=L"E_{%$(Ns[2]-1)}=%$(round2(Es[Ns[2]]))"
    ,xlabelsize=25,ylabelsize=25,titlesize=25,aspect=1)
ax3 = Axis(fig[2,1],xlabel=L"q",ylabel=L"p",title=L"E_{%$(Ns[3]-1)}=%$(round2(Es[Ns[3]]))"
    ,xlabelsize=25,ylabelsize=25,titlesize=25,aspect=1)
ax4 = Axis(fig[2,2],xlabel=L"q",ylabel=L"p",title=L"E_{%$(Ns[4]-1)}=%$(round2(Es[Ns[4]]))"
    ,xlabelsize=25,ylabelsize=25,titlesize=25,aspect=1)
h1 = heatmap!(ax1,qs,ps, Ws[:,:,Ns[1]], colormap = cgrad([:red,:white,:blue]),colorrange =(-lim,lim))
heatmap!(ax2,qs,ps, Ws[:,:,Ns[2]], colormap = cgrad([:red,:white,:blue]),colorrange =(-lim,lim))
heatmap!(ax3,qs,ps, Ws[:,:,Ns[3]], colormap = cgrad([:red,:white,:blue]),colorrange =(-lim,lim))
heatmap!(ax4,qs,ps, Ws[:,:,Ns[4]], colormap = cgrad([:red,:white,:blue]),colorrange =(-lim,lim))
cbar = Colorbar(fig[:,3], h1, width = 20, ticklabelsize = 20)
contour!(ax1,qs,ps, E, color = :black, levels = [Es[Ns[1]]], linewidth = 2.5)
contour!(ax2,qs,ps, E, color = :black, levels = [Es[Ns[2]]], linewidth = 2.5)
contour!(ax3,qs,ps, E, color = :black, levels = [Es[Ns[3]]], linewidth = 2.5)
contour!(ax4,qs,ps, E, color = :black, levels = [Es[Ns[4]]], linewidth = 2.5)
fig
##
save("waveFunc2WigFunc/Plots/harmonic_oscilator.png",fig)