using CairoMakie
##
include("fft_utils.jl")
include("solve1DSchrodinger.jl")
include("waveFunc2WigFunc.jl")
##
using SpecialFunctions,ClassicalOrthogonalPolynomials
ξ(χ,q) = exp(-q)/χ
s(χ,n) = 1/χ-2n-1
normalization(χ,n) = s(χ,n)*gamma(n+1)/gamma(1/χ-n)
N_max(χ) = Int(floor(.5*(1/χ-1)))
E(n,χ) = (n+.5)-χ*(n+.5)^2

function ψ(y,χ,n) 
	result = √(  normalization(χ,n)*ξ(χ,y)^s(χ,n)*exp(-ξ(χ,y))  )*laguerrel(n,s(χ,n),ξ(χ,y))
	ifelse(isnan(result),0.,result)
end
##
qs = LinRange(-log(2)-.1,8,1024)
ps_grid = DFTGrid(1024,10,1)
ps = direct_grid(ps_grid)
χ = .07
N_max(χ)
##
Ws = calculate_many_wigners((q,n)->ψ(q,χ,n),qs,ps_grid,7)
##
round2(x) = round(x,digits=2)
energies = [ χ*p^2+(1-exp(-q))^2/(4χ) for q in qs, p in ps ]
Ns = [1,3,5,7]
#lim = .8*maximum(Ws[:,:,n])
lim = .12

fig =Figure(fontsize=20,resolution=(1200,1000))
ax1 = Axis(fig[1,1],xlabel=L"q",ylabel=L"p",title=L"E_{%$(Ns[1]-1)}=%$( round2( E(Ns[1],χ) ) )"
    ,xlabelsize=25,ylabelsize=25,titlesize=25)
ax2 = Axis(fig[1,2],xlabel=L"q",ylabel=L"p",title=L"E_{%$(Ns[2]-1)}=%$(round2(E(Ns[2],χ)))"
    ,xlabelsize=25,ylabelsize=25,titlesize=25)
ax3 = Axis(fig[2,1],xlabel=L"q",ylabel=L"p",title=L"E_{%$(Ns[3]-1)}=%$(round2(E(Ns[3],χ)))"
    ,xlabelsize=25,ylabelsize=25,titlesize=25)
ax4 = Axis(fig[2,2],xlabel=L"q",ylabel=L"p",title=L"E_{%$(Ns[4]-1)}=%$(round2(E(Ns[4],χ)))"
    ,xlabelsize=25,ylabelsize=25,titlesize=25)
h1 = heatmap!(ax1,qs,ps, Ws[:,:,Ns[1]], colormap = cgrad([:red,:white,:blue]),colorrange =(-lim,lim))
heatmap!(ax2,qs,ps, Ws[:,:,Ns[2]], colormap = cgrad([:red,:white,:blue]),colorrange =(-lim,lim))
heatmap!(ax3,qs,ps, Ws[:,:,Ns[3]], colormap = cgrad([:red,:white,:blue]),colorrange =(-lim,lim))
heatmap!(ax4,qs,ps, Ws[:,:,Ns[4]], colormap = cgrad([:red,:white,:blue]),colorrange =(-lim,lim))
cbar = Colorbar(fig[:,3], h1, width = 20, ticklabelsize = 20)
contour!(ax1,qs,ps, energies, color = :black, levels = [ E(Ns[1],χ) ], linewidth = 2.5)
contour!(ax2,qs,ps, energies, color = :black, levels = [E(Ns[2],χ)], linewidth = 2.5)
contour!(ax3,qs,ps, energies, color = :black, levels = [E(Ns[3],χ)], linewidth = 2.5)
contour!(ax4,qs,ps, energies, color = :black, levels = [E(Ns[4],χ)], linewidth = 2.5)
fig
##
save("waveFunc2WigFunc/Plots/morse.png",fig)