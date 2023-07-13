includet("../../plot_config.jl")

V(x,y,μ) = (x^2/2-y)^2 + μ*x^2

V(x,y) = V(x,y,0.5)

xs = LinRange(-5,5,512)
ys = LinRange(-4,14,512)
Vs = [V(x,y) for y ∈ ys, x ∈ xs]
##
fig,ax,hm = heatmap(xs,ys,Vs,colorrange=(0,10))
Colorbar(fig[1, 2], hm)
fig