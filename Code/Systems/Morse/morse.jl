using SCCanonicalEnsemble
using FastGaussQuadrature
using SpecialFunctions
using ForwardDiff: derivative
using JLD2
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
χ = .12
θ_min = .01
θ_max = 7
N = 32
scatter_θs = LinRange(θ_min,θ_max,N)
line_θs = LinRange(θ_min,θ_max,4N)

Xs,ws = getNodesAndWeights(χ)
Us_ex = exact_energy(line_θs,χ,Es)
Us_sc = energy_quadrature(scatter_θs, χ, f!, H, Xs, ws)
Us_c = classical_energy.(line_θs,χ)
#jldsave("Results/Morse/energy_χ=$χ.jld2";scatter_θs,line_θs,χ,Us_ex,Us_sc,Us_c)


with_theme(theme) do
    f = Figure(;fonts = (; regular = texfont()))
    ax = Axis(f[1, 1],
    xlabel=L"\theta",xticks=0:5,
    ylabel = L"E")

    lines!(line_θs, Us_ex, color=colors[1],label="Quantum")
    scatter!(scatter_θs, first.(Us_sc), color=colors[2],label="Semiclassical",marker=:diamond)
    lines!(line_θs, Us_c, color=colors[3],label="Classical",linestyle=:dash)

    text!(ax, .42, .84, text=L"\chi = %$(try_conversion(χ))", space = :relative, fontsize=36)

    axislegend(; merge=true, position=:rt)
    f
    #save("Plots/Nelson/energy_μ=$μ.svg",f)
end
##
χ = .12
θ_min = .01
θ_max = 5
N = 32
scatter_θs = LinRange(θ_min,θ_max,N)
line_θs = LinRange(θ_min,θ_max,4N)

Cs_ex = exact_heat(line_θs,χ,Es)
Cs_sc = solve_equations(scatter_θs,χ,f!,χ-> getNodesAndWeights(χ,400),H,
output_func=heat_output,reduction=heat_reduction,callback=disc_caustic_callback)
Cs_c = classical_heat.(line_θs,χ)
##
p = plot(line_θs,Us_ex,
        ylabel=L"c/k",
        annotations = ((.9,.6), Plots.text(L"\chi=%$χ",10)),
        )   

    scatter!(scatter_θs,Us_sc,label=false,markershape=:diamond)
    plot!(line_θs,Us_c,label=false,line=:dot)
#jldsave("Results/Morse/heat_χ=$χ.jld2";scatter_θs,line_θs,χ,Cs_ex,Cs_sc,Cs_c)
##
using Plots,LaTeXStrings

default()
default(label=false,width=4,size=(354,250), markersize = 5, msw=0, 
palette=:Set1_3, tickfontsize=12,
guidefontsize=14, fontfamily="Computer Modern",dpi=1000,grid=false,framestyle = :box)

function make_plot(line_θs,scatter_θs,exact,sc,classical,χ,show_tick_legend)
    p = plot(line_θs,exact,
        ylabel=L"E/\hbar \omega",
        annotations = ((.8,.85), Plots.text(L"\chi=%$χ",14)),
        ylims=(0.01,min(first(exact)*1.05,2.4)))
    
    plot!(p,xlabel=L"\omega_0 \theta")

    scatter!(scatter_θs[1:2:end],sc[1:2:end],label=false,markershape=:diamond)
    plot!(line_θs,classical,label=false,line=:dot)
end

ps = []

for χ ∈ [.01,.04,.08,.12]
    _load(name) = load("Results/Morse/energy_χ=$χ.jld2",name)
    push!(ps,make_plot(_load("line_θs"),_load("scatter_θs"),
    _load("Us_ex"),_load("Us_sc"),_load("Us_c"),χ,length(ps)==3))
end

plot(ps[1],ps[2],ps[3],ps[4],size=(354*2,480),layout=(2,2),left_margin=3Plots.mm,)
#png("Plots/Morse/morse_energies")
##
using Plots,LaTeXStrings

default()
default(label=false,width=4,size=(354,250), markersize = 5, msw=0, 
palette=:Set1_3, tickfontsize=12,
guidefontsize=14, fontfamily="Computer Modern",dpi=1000,grid=false,framestyle = :box)

function make_plot(line_θs,scatter_θs,exact,sc,classical,χ,show_tick_legend)
    p = plot(line_θs,exact,
        ylabel=L"c/k",
        annotations = ((.45,.2), Plots.text(L"\chi=%$χ",14)),
        )
    
        plot!(p,xlabel=L"\omega_0 \theta")   

    scatter!(scatter_θs[1:2:end],sc[1:2:end],label=false,markershape=:diamond)
    plot!(line_θs,classical,label=false,line=:dot)
end

ps = []

for χ ∈ [.01,.04,.08,.12]
    _load(name) = load("Results/Morse/heat_χ=$χ.jld2",name)
    push!(ps,make_plot(_load("line_θs"),_load("scatter_θs"),
    _load("Cs_ex"),_load("Cs_sc"),_load("Cs_c"),χ,length(ps)==3))
end

plot(ps[1],ps[2],ps[3],ps[4],size=(354*2,480),layout=(2,2),left_margin=3Plots.mm)
#png("Plots/Morse/morse_heats")