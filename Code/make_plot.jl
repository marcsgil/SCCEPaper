include("plot_config.jl")
try_conversion(x) = isinteger(x) ? Int(x) : x
using JLD2
##
function make_plot(ylabel,pars,makepath, maketext, textpos)
    f = Figure(resolution=(1800, 1200), fontsize=48)
    positions = Iterators.product(1:2, 1:2) |> collect |> permutedims
    axs = [Axis(f[pos...];
        xticks=0:2:12, 
        ylabelsize=48, 
        xlabelsize=48, 
        xticklabelsize=40, 
        yticklabelsize=40) for pos ∈ positions]

    for n ∈ 1:4
        par = pars[n]
        path = makepath(par)
        line_θs = load(path, "line_θs")
        scatter_θs = load(path, "scatter_θs")
        quantum = load(path, "quantum")
        sc = load(path, "sc")
        classical = load(path, "classical")

        lines!(axs[n], line_θs, quantum, color=colors[1], linewidth=7)
        lines!(axs[n], scatter_θs, sc, color=colors[2], linestyle=:dot, linewidth=14)
        lines!(axs[n], line_θs, classical, color=colors[3], linewidth=7, linestyle=:dash)
        text!(axs[n], textpos..., text=maketext(par), space=:relative)
    end

    hidexdecorations!(axs[1, 1], ticks=false)
    hidexdecorations!(axs[2, 1], ticks=false)
    axs[1, 2].xlabel = L"\omega_0 \theta"
    axs[2, 2].xlabel = L"\omega_0 \theta"

    axs[1, 1].ylabel = ylabel
    axs[1, 2].ylabel = ylabel

    f,axs
end

ylabel = L"c / k"
pars = (0.1,0.3,0.6,1.)
makepath(χ) = "Results/Kerr/heat_χ=$χ.jld2"
maketext(χ) = L"\chi = %$(try_conversion(χ))"

f,axs = make_plot(ylabel,pars,makepath, maketext, (0.75,0.6))
ylims!(axs[1],0,.99)
ylims!(axs[2],0,.99)
ylims!(axs[3],0,.9)
ylims!(axs[4],0,.9)
f
#save("Plots/Kerr/heat_morse.pdf",f)