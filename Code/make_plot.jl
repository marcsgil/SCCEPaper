include("plot_config.jl")
using JLD2
##
function make_plot_2by2(ylabel,pars,makepath, maketext, textpos)
    f = Figure(resolution=(1800, 1200), fontsize=48)
    positions = Iterators.product(1:2, 1:2) |> collect |> permutedims
    axs = [Axis(f[pos...];
        xticks=0:7, ylabelsize=48, xlabelsize=48, xticklabelsize=40, yticklabelsize=40) for pos ∈ positions]

    for n ∈ 1:4
        par = pars[n]
        path = makepath(par)
        line_θs = load(path, "line_θs")
        scatter_θs = load(path, "scatter_θs")
        quantum = load(path, "quantum")
        sc = load(path, "sc")
        classical = load(path, "classical")

        lines!(axs[n], line_θs, quantum, color=colors[1], linewidth=7)
        scatter!(axs[n], scatter_θs, sc, color=colors[2], markersize=28)
        lines!(axs[n], line_θs, classical, color=colors[3], linewidth=7, linestyle=:dash)
        text!(axs[n], textpos..., text=maketext(par), space=:relative)
    end

    hidexdecorations!(axs[1, 1], ticks=false)
    hidexdecorations!(axs[2, 1], ticks=false)
    axs[1, 2].xlabel = L"\theta"
    axs[2, 2].xlabel = L"\theta"

    axs[1, 1].ylabel = ylabel
    axs[1, 2].ylabel = ylabel

    f,axs
end

function make_plot_2by1(ylabel,pars,makepath, maketext, textpos)
    f = Figure(resolution=(1800, 1200), fontsize=48)
    positions = Iterators.product(1:2, 1:1) |> collect 
    axs = [Axis(f[pos...];
        xticks=0:7, ylabelsize=48, xlabelsize=48, xticklabelsize=40, yticklabelsize=40) for pos ∈ positions]

    for n ∈ 1:3
        par = pars[n]
        path = makepath(par)
        line_θs = load(path, "line_θs")
        scatter_θs = load(path, "scatter_θs")
        quantum = load(path, "quantum")
        sc = load(path, "sc")
        classical = load(path, "classical")

        lines!(axs[n], line_θs, quantum, color=colors[1], linewidth=7)
        scatter!(axs[n], scatter_θs, sc, color=colors[2], markersize=28)
        lines!(axs[n], line_θs, classical, color=colors[3], linewidth=7, linestyle=:dash)
        text!(axs[n], textpos..., text=maketext(par), space=:relative)
    end

    hidexdecorations!(axs[1, 1], ticks=false)
    hidexdecorations!(axs[2, 1], ticks=false)
    axs[1, 2].xlabel = L"\theta"
    axs[2, 2].xlabel = L"\theta"

    axs[1, 1].ylabel = ylabel
    axs[1, 2].ylabel = ylabel

    f,axs
end

ylabel = L"E/\hbar \omega"
pars = (0.01,0.04,0.08,0.12)
makepath(χ) = "Results/Morse/energy_χ=$χ.jld2"
maketext(χ) = L"\chi = %$χ"

f,axs = make_plot(ylabel,pars,makepath, maketext, (0.75,0.6))
ylims!(axs[1],0,2.2)
ylims!(axs[2],0,2.2)
f
#save("Plots/Morse/energy_morse.pdf",f)