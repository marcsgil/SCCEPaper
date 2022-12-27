using GLMakie

V(r,λ=1) = (r[1]^2+r[2]^2)/2 + λ*(r[1]^4*r[2]^2)
##
rmax = 6
N=128
rs = LinRange(-rmax,rmax,N)
zs = map(V,Iterators.product(rs,rs))
##
fig = Figure()
ax = Axis3(fig[1,1])
#zlims!(ax,-5,5)
GLMakie.surface!(ax,rs,rs,zs)
fig