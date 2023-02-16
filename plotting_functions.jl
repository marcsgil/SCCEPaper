using CairoMakie,MakiePublication

set_theme!(theme_aps())

function comparison_plot(θs_lines,θs_scatter,Us_ex,Us_sc;title="")
    CairoMakie.activate!()
    f = CairoMakie.Figure(fontsize=24)
    ax = CairoMakie.Axis(f[1, 1],
        xlabel = L"θ",
        ylabel = "Energy",
        title = title
    )
    ylims!(ax,last(Us_ex)*2^sign(-last(Us_ex)),first(Us_ex)*1.05,)
    lines!(ax,θs_lines,Us_ex,label="Exact",color=:black)
    scatter!(ax,θs_scatter,Us_sc,label="SC",color=:red)
    axislegend()
    f
end

function comparison_plot(θs_lines,θs_scatter,Us_ex,Us_sc,Us_c;title="")
    CairoMakie.activate!()
    f = CairoMakie.Figure(fontsize=24)
    ax = CairoMakie.Axis(f[1, 1],
        xlabel = L"θ",
        ylabel = "Energy",
        title = title,
    )
    ylims!(ax,last(Us_ex)*2^sign(-last(Us_ex)),first(Us_ex)*1.05,)
    lines!(ax,θs_lines,Us_ex,label="Exact")
    lines!(ax,θs_scatter,Us_sc,label="Semi-classical")
    lines!(ax,θs_lines,Us_c,label="Classical")
    axislegend()
    f
end

#=using GLMakie
function analysis_plot(sols,ps,qs,θs,χ)
    GLMakie.activate!()

    fig = Figure(resolution=(1600,700), fontsize=24)
    ax1 = GLMakie.Axis(fig[1,1],xlabel="q",ylabel="p",title="Integrando F. de Partição")
    ax2 = GLMakie.Axis(fig[1,3],xlabel="q",ylabel="p",title="Determinante Jacobiano")

    sl = Slider(fig[2, :], range = eachindex(θs), startvalue = 1)

    title = Label(fig[0, :], L"\theta = %$(round(θs[1],digits=2)),\chi=%$χ", fontsize = 30)
    lift(sl.value) do n
        title.text = L"\theta = %$(round(θs[n],digits=2)),\chi=%$χ"
    end

    N = length(ps)

    Zs = lift(n->reshape(view(sols,1,n,:),N,N) |> transpose,sl.value)
    dets = lift(n->reshape(view(sols,2,n,:),N,N) |> transpose,sl.value)

    hm1 = heatmap!(ax1,qs,ps,Zs,colormap=:afmhot,colorrange=(0,3))
    Colorbar(fig[1, 2], hm1)

    hm2 = heatmap!(ax2,qs,ps,dets,colormap=:redsblues,colorrange=(-1,1))
    Colorbar(fig[1, 4], hm2)

    fig
end=#