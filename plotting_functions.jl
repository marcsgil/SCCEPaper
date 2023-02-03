using CairoMakie

function comparison_plot(θs_ex,θs_sc,Us_ex,Us_sc,title="")
    CairoMakie.activate!()
    f = CairoMakie.Figure(fontsize=24)
    ax = CairoMakie.Axis(f[1, 1],
        xlabel = L"θ",
        xlabelsize=32,
        ylabel = "Energy",
        ylabelsize=24,
        title = title,
        titlesize = 32
    )
    lines!(ax,θs_ex,Us_ex,label="Exact",color=:black)
    scatter!(ax,θs_sc,Us_sc,label="SC",color=:red,marker=:diamond)
    axislegend()
    f
end

#=

##
using GLMakie
GLMakie.activate!()
##
N = 64
χ=.12
ps = LinRange(-1/2χ,1/2χ,N)
qs = LinRange(-log(2),1,N)
xs = [SA[p,q] for p in ps, q in qs]

θs = LinRange(0,3,64)


sols = solve_equations(θs,χ,fy,fx,(par)->(xs,ones(N^2)),H,
abstol=1e-12,reltol=1e-12,alg=Vern8(),output_func=analysis_output,callback=strong_callback)
##
fig = Figure(resolution=(1600,700), fontsize=24)
ax1 = GLMakie.Axis(fig[1,1],xlabel="q",ylabel="p",title="Integrando F. de Partição",xticks=-5:1:5)
ax2 = GLMakie.Axis(fig[1,3],xlabel="q",ylabel="p",title="Determinante Jacobiano",xticks=-5:1:5)

sl = Slider(fig[2, :], range = eachindex(θs), startvalue = 1)

title = Label(fig[0, :], L"\theta = %$(round(θs[1],digits=2)),\chi=%$χ", fontsize = 30)
lift(sl.value) do n
    title.text = L"\theta = %$(round(θs[n],digits=2)),\chi=%$χ"
end

Zs = lift(n->reshape(view(sols,1,n,:),N,N) |> transpose,sl.value)
dets = lift(n->reshape(view(sols,2,n,:),N,N) |> transpose,sl.value)

hm1 = heatmap!(ax1,qs,ps,Zs,colormap=:afmhot,colorrange=(0,3))
Colorbar(fig[1, 2], hm1)

hm2 = heatmap!(ax2,qs,ps,dets,colormap=:redsblues,colorrange=(-1,1))
Colorbar(fig[1, 4], hm2)
=#