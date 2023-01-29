using CairoMakie

function comparison_plot(θs_ex,θs_sc,Us_ex,Us_sc,title)
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