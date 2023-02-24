using ExactCanonicalEnsemble,SCCanonicalEnsemble
using ClassicalCanonicalEnsemble
using SolveSchrodinger
using ProgressMeter,JLD2
##
V(q,μ) = (q[1]^2/2-q[2])^2 + μ*q[1]^2
H(x,μ) = (x[1]^2+x[2]^2)/2 + (x[3]^2/2-x[4])^2 + μ*x[3]^2
f! = get_equations_of_motion(H,2)
##
μ = .5

Es,ψs = solveSchrodinger(xs,ys,V;nev=60,par=μ)
ex_U(θ,Es) = sum(E->E*exp(-θ*E),Es)/sum(E->exp(-θ*E),Es)

θ = 1
ex_U(θ,Es)

energyMonteCarlo(θ,μ,H,2,f!,10^5,callback=disc_caustic_callback,alg=Tsit5(),reltol=1e-3,abstol=1e-6)
##
μ = 2
θ_min = .1
θ_max = 3
N = 16
scatter_θs = LinRange(θ_min,θ_max,N)
line_θs = LinRange(θ_min,θ_max,4N)

xs = LinRange(-4.5,4.5,160)
ys = LinRange(-4,5,160)
Es,ψs = solveSchrodinger(xs,ys,V;nev=70,par=μ)
Us_ex = [exact_energy(θ,Es) for θ in line_θs]
Us_sc = @showprogress [energyMonteCarlo(θ,μ,H,2,f!,10^6,callback=disc_caustic_callback,alg=Tsit5(),reltol=1e-3,abstol=1e-6) for θ in scatter_θs]
Us_c = classical_energy(line_θs,2,μ,potential=V,alg=CubatureJLh(),reltol=1e-3)
jldsave("Results/Nelson/energy_μ=$μ.jld2";scatter_θs,line_θs,μ,Us_ex,Us_sc,Us_c)
##
μ = 2
θ_min = .1
θ_max = 3
N = 16
scatter_θs = LinRange(θ_min,θ_max,N)
line_θs = LinRange(θ_min,θ_max,4N)

xs = LinRange(-4.5,4.5,160)
ys = LinRange(-4,5,160)
Es,ψs = solveSchrodinger(xs,ys,V;nev=70,par=μ)
Us_ex = [exact_heat(θ,Es) for θ in line_θs]
Us_sc = @showprogress [heatMonteCarlo(θ,μ,H,2,f!,10^6,callback=disc_caustic_callback,alg=Tsit5(),reltol=1e-3,abstol=1e-6) for θ in scatter_θs]
Us_c = classical_energy(line_θs,2,μ,potential=V,alg=CubatureJLh(),reltol=1e-3)
jldsave("Results/Nelson/energy_μ=$μ.jld2";scatter_θs,line_θs,μ,Us_ex,Us_sc,Us_c)
##
using Plots,LaTeXStrings

default(label=false,width=3,size=(354,250), markersize = 3, msw=0, 
palette=:Set1_3, tickfontsize=8,
guidefontsize=8, fontfamily="Computer Modern",dpi=1000)

p = plot(line_θs,Us_ex,
        ylabel=L"U/\hbar \omega",xlabel=L"\omega \theta",
        annotations = ((.9,.9), Plots.text(L"\mu=%$μ",10)),
        yaxis=:log) 

plot!(line_θs,Us_c,label=false,line=:dot)
scatter!(scatter_θs,Us_sc,label=false,markershape=:diamond)

##
#png("Plots/Nelson/μ=$μ")
##