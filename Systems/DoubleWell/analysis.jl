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
N = 256
ps = LinRange(-12,12,N)
qs = LinRange(-6,6,N)

θs = LinRange(0,2,120)
χ=2
##
sols = solve_equations(θs,χ,fy,fx,(par)->([SA[p,q] for p in ps, q in qs],ones(N^2)),
abstol=1e-12,reltol=1e-12,alg=Vern8(),output_func=analysis_output,stop_at_caustic=false)
##
fig = Figure(resolution=(1600,700), fontsize=24)
ax1 = GLMakie.Axis(fig[1,1],xlabel="q",ylabel="p",title="Integrando F. de Partição")
ax2 = GLMakie.Axis(fig[1,3],xlabel="q",ylabel="p",title="Determinante Jacobiano")

sl = Slider(fig[2, :], range = eachindex(θs), startvalue = 1)

title = Label(fig[0, :], L"\theta = %$(round(θs[1],digits=2)),\chi=%$χ", textsize = 30)
lift(sl.value) do n
    title.text = L"\theta = %$(round(θs[n],digits=2)),\chi=%$χ"
end

Zs = lift(n->reshape(view(sols,1,n,:),N,N) |> transpose,sl.value)
dets = lift(n->reshape(view(sols,2,n,:),N,N) |> transpose,sl.value)

hm1 = heatmap!(ax1,qs,ps,Zs,colormap=:afmhot,colorrange=(0,3))
Colorbar(fig[1, 2], hm1)

hm2 = heatmap!(ax2,qs,ps,dets,colormap=:redsblues,colorrange=(-1,1))
Colorbar(fig[1, 4], hm2)

fig
##
fig = Figure(resolution=(1600,700), fontsize=24)
ax1 = GLMakie.Axis(fig[1,1],xlabel="q",ylabel="p",title="Integrando F. de Partição")
ax2 = GLMakie.Axis(fig[1,3],xlabel="q",ylabel="p",title="Determinante Jacobiano")
title = Label(fig[0, :], L"\theta = %$(round(θs[1],digits=2)),\chi=%$χ", textsize = 30)

function update!(fig,title,n)
    title.text = L"\theta = %$(round(θs[n],digits=2)),\chi=%$χ"

    hm1 = heatmap!(ax1,qs,ps,reshape(view(sols,1,n,:),N,N) |> transpose,colormap=:afmhot,colorrange=(0,2))
    n == 1 ? Colorbar(fig[1, 2], hm1) : nothing

    hm2 = heatmap!(ax2,qs,ps,reshape(view(sols,2,n,:),N,N) |> transpose,colormap=:redsblues,colorrange=(-1,1))
    n == 1 ? Colorbar(fig[1, 4], hm2) : nothing
end

record(fig, "Plots/double_well.mp4",framerate=18) do io
    for n in eachindex(θs)
        update!(fig,title,n)     # animate figure
        recordframe!(io)  # record a new frame
    end
end