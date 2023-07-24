includet("../../plot_config.jl")

V(x,y,μ) = (x^2/2-y)^2 + μ*x^2

V(x,y) = V(x,y,0.5)

xs = LinRange(-2,2,512)
ys = LinRange(-2,2,512)
Vs = [V(x,y) for y ∈ ys, x ∈ xs]
##
fig,ax,hm = surface(xs,ys,Vs,colorrange=(0,1))
zlims!(ax,0,1)
fig