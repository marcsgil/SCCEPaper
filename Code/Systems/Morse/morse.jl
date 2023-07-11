using SCCanonicalEnsemble, QuantumCanonicalEnsemble
using FastGaussQuadrature
using SpecialFunctions
using ForwardDiff: derivative
using JLD2
includet("../../plot_config.jl")
##
H(x,χ) = χ*x[1]^2+(1-exp(-x[2]))^2/(4χ)
H2(x,χ) = H(x,χ)^2 - (2exp(-2x[2])-exp(-x[2]))/4
Es(χ) = [ (n+0.5)-χ*(n+0.5)^2 for n in 0:floor(Int64, 0.5*(1/χ-1))]

function f!(du,u,χ)
    du[1] = -4χ*u[3]
    du[2] = ( -exp(-u[4])*cos(0.5*u[1]) + exp(-2*u[4])*cos(u[1]) )/χ
    du[3] = (exp(-u[4])*sin(0.5*u[1]) - exp(-2*u[4])*sin(u[1]) )/(2χ)
    du[4] = -χ*u[2]
end

function getNodesAndWeights(χ,N=300)
    Ps,wPs = gausslegendre(N)
    Qs,wQs = gausschebyshev(N,3)

    [[√(1-Q^2)*P/(2χ),-log(1-Q)] for P in Ps, Q in Qs] |> vec |> stack, [w1*w2 for w1 in wPs, w2 in wQs] |> vec
end

Z(θ,χ) = dawson(√(θ/4χ))/√θ

function classical_energy(θ,par)
    -derivative( θ->log(Z(θ,par)), θ )
end

function classical_heat(θ,par)
    θ^2*derivative(θ->derivative( θ->log(Z(θ,par)), θ ),θ)
end
##
χ = .01
θ_min = .4
θ_max = 7
N = 16
scatter_θs = LinRange(θ_min,θ_max,N)
line_θs = LinRange(θ_min,θ_max,4N)

Xs,ws = getNodesAndWeights(χ)
Us_qu = quantum_energy(line_θs,χ,Es)
Us_sc = energy_quadrature(scatter_θs, χ, f!, H, Xs, ws)
Us_c = classical_energy.(line_θs,χ)

with_theme(theme) do
    f = Figure(;fonts = (; regular = texfont()))
    ax = Axis(f[1, 1],
    xlabel=L"\omega \theta",xticks=0:θ_max,
    ylabel = L"E / \hbar \omega")

    lines!(line_θs, Us_qu, color=colors[1],label="Quantum")
    scatter!(scatter_θs, Us_sc, color=colors[2],label="Semiclassical",marker=:diamond)
    lines!(line_θs, Us_c, color=colors[3],label="Classical",linestyle=:dash)

    text!(ax, .42, .84, text=L"\chi = %$χ", space = :relative, fontsize=36)

    axislegend(; merge=true, position=:rt)
    f
    #save("Plots/Morse/energy_χ=$χ.svg",f)
    #jldsave("Results/Morse/energy_χ=$χ.jld2";scatter_θs,line_θs,χ,Us_qu,Us_sc,Us_c)
end
##
for χ ∈ (.01,.04,.08,.12)
    θ_min = .01
    θ_max = 5
    N = 16
    scatter_θs = LinRange(θ_min,θ_max,N)
    line_θs = LinRange(θ_min,θ_max,8N)
    
    Xs,ws = getNodesAndWeights(χ)
    Cs_qu = quantum_heat(line_θs,χ,Es)
    Cs_sc = heat_quadrature(scatter_θs, χ, f!, H, Xs, ws)
    Cs_c = classical_heat.(line_θs,χ)
    jldsave("Results/Morse/heat_χ=$χ.jld2";scatter_θs,line_θs,χ,Us_qu,Us_sc,Us_c)
end

with_theme(theme) do
    f = Figure(;fonts = (; regular = texfont()))
    ax = Axis(f[1, 1],
    xlabel=L"\omega \theta",xticks=0:θ_max,
    ylabel = L"c / k")

    lines!(line_θs, Cs_qu, color=colors[1],label="Quantum")
    scatter!(scatter_θs, first.(Cs_sc), color=colors[2],label="Semiclassical",marker=:diamond)
    lines!(line_θs, Cs_c, color=colors[3],label="Classical",linestyle=:dash)

    text!(ax, .75, .85, text=L"\chi = %$χ", space = :relative, fontsize=36)

    axislegend(; merge=true, position=(1,.6))
    f
    #save("Plots/Morse/energy_χ=$χ.svg",f)
    #jldsave("Results/Morse/energy_χ=$χ.jld2";scatter_θs,line_θs,χ,Us_qu,Us_sc,Us_c)
end