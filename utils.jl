using CairoMakie

function make_plot(xline,yline,xscatter,yscatter,anotation,ylabel)
    f = Figure(fontsize=24)
    ax = Axis(f[1, 1])
    lines!(ax,xline,yline,label="Exact",color=:black)
    scatter!(ax,xscatter,yscatter,label="SC",color=:red,marker=:diamond)
    axislegend()
    px = .45
    py = .91
    text!(ax,px*xline[end]+(1-px)xline[1],py*yscatter[1]+(1-py)yscatter[end],text=anotation,textsize=36)
    f
end