using StaticArrays,GLMakie

include("../../double_phase_space.jl")
include("double_well.jl")

function analysis_output(sol,i,θs,par,node,weight)
    output = Array{eltype(θs)}(undef,2,length(θs))
    E = H(node,par)
    for (n,u) in enumerate(sol.u)
        output[1,n] = exp( u.Δ - θs[n]*E )*√abs(det(u.jac_x))
        output[2,n] = det(u.jac_x)
    end
    output,false
end
##
N = 128
ps = LinRange(-12,12,N)
qs = LinRange(-6,6,N)

θs = LinRange(0,3,64)
χ=4
##
sols = solve_equations(θs,χ,fy,fx,(par)->([SA[p,q] for p in ps, q in qs],ones(N^2)),
abstol=1e-12,reltol=1e-12,alg=Vern8(),stop_at_caustic=false,output_func=analysis_output)
##
fig = Figure(resolution=(1600,600))
ax1 = GLMakie.Axis(fig[1,1])
ax2 = GLMakie.Axis(fig[1,3])

sl = Slider(fig[2, :], range = eachindex(θs), startvalue = 1)

title = Label(fig[0, :], L"\theta = %$(round(θs[1],digits=2))", textsize = 30)
lift(sl.value) do n
    title.text = L"\theta = %$(round(θs[n],digits=2))"
end

Zs = lift(n->reshape(view(sols,1,n,:),N,N) |> transpose,sl.value)
dets = lift(n->reshape(view(sols,2,n,:),N,N) |> transpose,sl.value)

hm1 = heatmap!(ax1,qs,ps,Zs,colormap=:hot,colorrange=(0,1.5))
Colorbar(fig[1, 2], hm1)

hm2 = heatmap!(ax2,qs,ps,dets,colormap=:redsblues,colorrange=(-1,1))
Colorbar(fig[1, 4], hm2)

fig