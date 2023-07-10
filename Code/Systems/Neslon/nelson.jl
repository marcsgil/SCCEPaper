using ExactCanonicalEnsemble, SCCanonicalEnsemble, ClassicalCanonicalEnsemble
using SolveSchrodinger
using ProgressMeter, CairoMakie, MathTeXEngine, ColorSchemes
##
V(q, μ) = (q[1]^2 / 2 - q[2])^2 + μ * q[1]^2
H(x, μ) = (x[1]^2 + x[2]^2) / 2 + (x[3]^2 / 2 - x[4])^2 + μ * x[3]^2
f! = get_equations_of_motion(H, 2)

##
μ = 1.
xs = LinRange(-4.5, 4.5, 160)
ys = LinRange(-4, 5, 160)
Es, ψs = solveSchrodinger(xs, ys, V; nev=70, par=μ)
θs = LinRange(0.1,5,24)
sc_quantities = Vector{NTuple{2,eltype(θs)}}(undef, length(θs))

@showprogress for (n, θ) ∈ enumerate(θs)
    lb = [30 * exp(-θ) for _ ∈ 1:4]
    ub = -lb
    sc_quantities[n] = heat_integrals(θ, μ, f!, H, ub, lb)
end

Us_sc = first.(sc_quantities)
Cs_sc = last.(sc_quantities)
Us_ex = exact_energy(θs, Es)
Cs_ex = exact_heat(θs,Es)
Us_c = [2/θ for θ ∈ θs]
##
colors = colorschemes[:seaborn_bright]
theme = Theme(
    fontsize=32,
    Axis=(xlabelsize=36, xlabelpadding=0,
    ylabelsize=36, ylabelpadding=0,
    xticklabelsize = 28, yticklabelsize=28,
    xgridvisible=false, ygridvisible=false,
    xtickalign=1, ytickalign=1,
    yticksize=12, xticksize=12),
    Lines = (linewidth = 5,),
    Scatter = (markersize=20,),
    Legend = (labelsize=28,)
)
try_conversion(x) = isinteger(x) ? Int(x) : x

with_theme(theme) do
    f = Figure(;fonts = (; regular = texfont()))
    ax = Axis(f[1, 1],
    xlabel=L"\theta",xticks=0:5,
    ylabel = L"E", yscale=log2,yticks =[2.0^n for n ∈ -1:4])

    lines!(θs, Us_ex, color=colors[1],label="Quantum")
    scatter!(θs, first.(Us_sc), color=colors[2],label="Semiclassical",marker=:diamond)
    lines!(θs, Us_c, color=colors[3],label="Classical",linestyle=:dash)

    text!(ax, .42, .84, text=L"\mu = %$(try_conversion(μ))", space = :relative, fontsize=36)

    axislegend(; merge=true, position=:rt)
    f
    #save("Plots/Nelson/energy_μ=$μ.svg",f)
end

with_theme(theme) do
    f = Figure(;fonts = (; regular = texfont()))
    ax = Axis(f[1, 1],
    xlabel=L"\theta",xticks=0:0.5:2, ylabel = L"c/k")

    lines!(θs, Cs_ex, color=colors[1],label="Quantum")
    scatter!(θs, first.(Cs_sc), color=colors[2],label="Semiclassical",marker=:diamond)
    hlines!(2, color=colors[3],label="Classical",linestyle=:dash,linewidth = 5)

    text!(ax, .8, .6, text=L"\mu = %$(try_conversion(μ))", space = :relative, fontsize=36)

    axislegend(; merge=true, position=(1,0.9))
    xlims!(ax,0,2.1)
    f
    save("Plots/Nelson/heat_μ=$μ.svg",f)
end