using CairoMakie

function comparison_plot(θs_ex,θs_sc,Us_ex,Us_sc)
    CairoMakie.activate!()
    f = CairoMakie.Figure(fontsize=24)
    ax = CairoMakie.Axis(f[1, 1],
        xlabel = L"θ",
        xlabelsize=32,
        ylabel = "Energy",
        ylabelsize=24
    )
    lines!(ax,θs_ex,Us_ex,label="Exact",color=:black)
    scatter!(ax,θs_sc,Us_sc,label="SC",color=:red,marker=:diamond)
    axislegend()
    px = .45
    py = .91
    text!(ax,px*θs_ex[end]+(1-px)θs_ex[1],py*Us_ex[1]+(1-py)Us_ex[end],text=L"\chi=%$χ",fontsize=36)
    f
end