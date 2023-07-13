using SCCanonicalEnsemble, QuantumCanonicalEnsemble, SolveSchrodinger
using JLD2, ProgressMeter
includet("../../plot_config.jl")
try_conversion(x) = isinteger(x) ? Int(x) : x
##
V(q, μ) = (q[1]^2 / 2 - q[2])^2 + μ * q[1]^2
H(x, μ) = (x[1]^2 + x[2]^2) / 2 + (x[3]^2 / 2 - x[4])^2 + μ * x[3]^2
f! = get_equations_of_motion(H, 2)
##
μ = .5
xs = LinRange(-4.5, 4.5, 160)
ys = LinRange(-4, 5, 160)
Es, ψs = solveSchrodinger(xs, ys, V; nev=70, par=μ)
θ_min = .2
θ_max = 5
N = 16
scatter_θs = LinRange(θ_min,θ_max,N)
line_θs = LinRange(θ_min,θ_max,4N)
Us_sc = Vector{eltype(scatter_θs)}(undef, length(scatter_θs))
##
θ = 1. 
ub = [30 * exp(-θ) for _ ∈ 1:4]
lb = -ub
quantum_energy(θ, Es)
energy_integrals(θ, μ, f!, H, lb, ub, reltol=1e-2,abstol=1e-3,maxiters=10^5)
##
ub = [5,5,5,5]
lb = [-5,-5,-5,-1]
energy_integrals2(θ, μ, f!, H, lb,ub, reltol=1e-2,abstol=1e-3,maxiters=10^5)
##
@showprogress for (n, θ) ∈ enumerate(scatter_θs)
    lb = [30 * exp(-θ) for _ ∈ 1:4]
    ub = -lb
    Us_sc[n] = energy_integrals(θ, μ, f!, H, ub, lb, reltol=1e-2,abstol=1e-3,maxiters=10^5)
end

Us_qu = quantum_energy(line_θs, Es)
Us_c = [2/θ for θ ∈ line_θs]
##

with_theme(theme) do
    f = Figure(;fonts = (; regular = texfont()))
    ax = Axis(f[1, 1],
    xlabel=L"\theta",xticks=0:5,
    ylabel = L"E", yscale=log2,yticks =[2.0^n for n ∈ -1:4])

    lines!(line_θs, Us_qu, color=colors[1],label="Quantum")
    scatter!(scatter_θs, first.(Us_sc), color=colors[2],label="Semiclassical",marker=:diamond)
    lines!(line_θs, Us_c, color=colors[3],label="Classical",linestyle=:dash)

    text!(ax, .42, .84, text=L"\mu = %$(try_conversion(μ))", space = :relative, fontsize=36)

    axislegend(; merge=true, position=:rt)
    f
    #jldsave("Results/Nelson/energy_μ=$μ.jld2";scatter_θs,line_θs,μ,Us_qu,Us_sc,Us_c)
end
##
θ_min = .2
θ_max = 2
N = 16
scatter_θs = LinRange(θ_min,θ_max,N)
line_θs = LinRange(.001,θ_max,4N)
Cs_sc = Vector{eltype(scatter_θs)}(undef, length(scatter_θs))

##
θ = scatter_θs[1]
ub = [12,12,60,12]
lb = [-12,-12,-10,-12]
heat_integrals(θ, μ, f!, H, ub, lb, reltol=1e-3,abstol=1e-3,maxiters=10^5)
##
data = heat_integrals(θ, μ, f!, H, ub, lb, reltol=1e-2,abstol=1e-3,maxiters=10^5)
f,ax,hm = heatmap(view(reshape(data,3,128,128),3,:,:))
Colorbar(f[1, 2], hm)
f
##
@showprogress for (n, θ) ∈ enumerate(scatter_θs)
    ub = [80 for _ ∈ 1:4] * exp(-θ)
    lb = -ub
    Cs_sc[n] = heat_integrals(θ, μ, f!, H, ub, lb, reltol=1e-2,abstol=1e-3,maxiters=10^5)
end

Cs_qu = quantum_heat(line_θs, Es)
Cs_c = [2 for θ ∈ line_θs]

with_theme(theme) do
    f = Figure(;fonts = (; regular = texfont()))
    ax = Axis(f[1, 1],
    xlabel=L"\theta", ylabel = L"c/k")

    lines!(line_θs, Cs_qu, color=colors[1],label="Quantum")
    scatter!(scatter_θs, first.(Cs_sc), color=colors[2],label="Semiclassical",marker=:diamond)
    lines!(line_θs,Cs_c, color=colors[3],label="Classical",linestyle=:dash)

    #text!(ax, .8, .6, text=L"\mu = %$(try_conversion(μ))", space = :relative, fontsize=36)

    #axislegend(; merge=true, position=(1,0.9))
    f
    #jldsave("Results/Nelson/heat_μ=$μ.jld2";scatter_θs,line_θs,μ,Us_qu,Us_sc,Us_c)
end
