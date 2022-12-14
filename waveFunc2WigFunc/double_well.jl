using CairoMakie
##
include("../fft_utils.jl")
include("../waveFunc2WigFunc.jl")
##
V(q,χ) = χ*q^4/4+(1/2-χ)*q^2
qs = LinRange(-6,6,1024)
χ = 1
Es,ψs = solveSchrodinger(qs,V,par=χ)
lines(qs,ψs[:, 21])
lines(qs,abs2.(ψs[:, 22]))
Es[1]
##
ps_grid = DFTGrid(1024,6,1)
ps = direct_grid(ps_grid)
##
Ws = calculate_many_wigners(ψs,qs,ps_grid,16)
##
round2(x) = round(x,digits=2)
E = [p^2/2 + V(q,χ) for q in qs, p in ps]
Ns = [1,2,3,4]
lim = .8*maximum(Ws[:,:,1])
#lim = .12
fig =Figure(fontsize=20,resolution=(1200,1000))
ax1 = CairoMakie.Axis(fig[1,1],xlabel=L"q",ylabel=L"p",title=L"E_{%$(Ns[1]-1)}=%$(round2(Es[Ns[1]]))"
    ,xlabelsize=25,ylabelsize=25,titlesize=25)
ax2 = CairoMakie.Axis(fig[1,2],xlabel=L"q",ylabel=L"p",title=L"E_{%$(Ns[2]-1)}=%$(round2(Es[Ns[2]]))"
    ,xlabelsize=25,ylabelsize=25,titlesize=25)
ax3 = CairoMakie.Axis(fig[2,1],xlabel=L"q",ylabel=L"p",title=L"E_{%$(Ns[3]-1)}=%$(round2(Es[Ns[3]]))"
    ,xlabelsize=25,ylabelsize=25,titlesize=25)
ax4 = CairoMakie.Axis(fig[2,2],xlabel=L"q",ylabel=L"p",title=L"E_{%$(Ns[4]-1)}=%$(round2(Es[Ns[4]]))"
    ,xlabelsize=25,ylabelsize=25,titlesize=25)
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
save("waveFunc2WigFunc/Plots/double_well.png",fig)