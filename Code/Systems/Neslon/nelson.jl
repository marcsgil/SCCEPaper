using SCCanonicalEnsemble, QuantumCanonicalEnsemble, SolveSchrodinger
using JLD2, ProgressMeter
includet("../../plot_config.jl")
##
V(q, μ) = (q[1]^2 / 2 - q[2])^2 + μ * q[1]^2
H(x, μ) = (x[1]^2 + x[2]^2) / 2 + (x[3]^2 / 2 - x[4])^2 + μ * x[3]^2
f! = get_equations_of_motion(H, 2)
##
μ = 1.0
xs = LinRange(-4.5, 4.5, 160)
ys = LinRange(-4, 5, 160)
Es, ψs = solveSchrodinger(xs, ys, V; nev=70, par=μ)
θ_min = .2
θ_max = 5
N = 16
scatter_θs = LinRange(θ_min,θ_max,N)
line_θs = LinRange(θ_min,θ_max,4N)
sc = Vector{eltype(scatter_θs)}(undef, length(scatter_θs))
##
θ = 1. 
ub = [45 * exp(-θ) for _ ∈ 1:4]
lb = -ub
quantum_energy(θ, Es)
energy_integrals(θ, μ, f!, H, lb, ub, reltol=1e-6,abstol=1e-6,maxiters=10^5, dif_eq_alg=Vern6())
#2/θ
##
@showprogress for (n, θ) ∈ enumerate(scatter_θs)
    lb = [30 * exp(-θ) for _ ∈ 1:4]
    ub = -lb
    sc[n] = energy_integrals(θ, μ, f!, H, ub, lb, reltol=1e-6,abstol=1e-6,maxiters=10^5, dif_eq_alg=Vern6())
end

quantum = quantum_energy(line_θs, Es)
classical = [2/θ for θ ∈ line_θs]
##

with_theme(theme) do
    f = Figure(;fonts = (; regular = texfont()))
    ax = Axis(f[1, 1],
    xlabel=L"\theta",xticks=0:5,
    ylabel = L"E", yscale=log2,yticks =[2.0^n for n ∈ -1:4])

    lines!(line_θs, quantum, color=colors[1],label="Quantum")
    scatter!(scatter_θs, first.(sc), color=colors[2],label="Semiclassical",marker=:diamond)
    lines!(line_θs, classical, color=colors[3],label="Classical",linestyle=:dash)

    text!(ax, .42, .84, text=L"\mu = %$(try_conversion(μ))", space = :relative, fontsize=36)

    axislegend(; merge=true, position=:rt)
    f
    jldsave("Results/Nelson/energy_μ=$μ.jld2";scatter_θs,line_θs,μ,quantum,sc,classical)
end
##
θ_min = .1
θ_max = 2
N = 16
scatter_θs = LinRange(θ_min,θ_max,N)
line_θs = LinRange(θ_min,θ_max,4N)
sc = Vector{eltype(scatter_θs)}(undef, length(scatter_θs))
##
θ = scatter_θs[1]
#ub = [4 , 4, 4, 4.] * exp(1.6*(2-θ))
#lb = [-4 , -4, -4, -1.5] * exp(1.6*(2-θ))
ub = [60 * exp(-θ) for _ ∈ 1:4]
lb = -ub
quantum_heat(θ, Es)
data = heat_integrals(θ, μ, f!, H, lb, ub, reltol=1e-8,abstol=1e-8,maxiters=10^6, dif_eq_alg=Vern8())

fig,ax,hm = heatmap(reshape(data[3,:],128,128))
maximum(data[3,:])
##
@showprogress for (n, θ) ∈ enumerate(scatter_θs)
    ub = [60 * exp(-θ) for _ ∈ 1:4]
    lb = -ub
    sc[n] = heat_integrals(θ, μ, f!, H, lb, ub, reltol=1e-8,abstol=1e-8,maxiters=10^6, dif_eq_alg=Vern8())
end

quantum = quantum_heat(line_θs, Es)
classical = [2 for θ ∈ line_θs]

with_theme(theme) do
    f = Figure(;fonts = (; regular = texfont()))
    ax = Axis(f[1, 1],
    xlabel=L"\theta", ylabel = L"c/k")

    lines!(line_θs, quantum, color=colors[1],label="Quantum")
    scatter!(scatter_θs, sc, color=colors[2],label="Semiclassical",marker=:diamond)
    lines!(line_θs,classical, color=colors[3],label="Classical",linestyle=:dash)

    #text!(ax, .8, .6, text=L"\mu = %$(try_conversion(μ))", space = :relative, fontsize=36)

    #axislegend(; merge=true, position=(1,0.9))
    f
    #jldsave("Results/Nelson/heat_μ=$μ.jld2";scatter_θs,line_θs,μ,quantum,sc,classical)
end
