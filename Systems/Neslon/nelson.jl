using ExactCanonicalEnsemble,SCCanonicalEnsemble
using ClassicalCanonicalEnsemble
using SolveSchrodinger
using ProgressMeter,JLD2,OnlineStats
using ThreadsX
##
V(q,μ) = (q[1]^2/2-q[2])^2 + μ*q[1]^2
H(x,μ) = (x[1]^2+x[2]^2)/2 + (x[3]^2/2-x[4])^2 + μ*x[3]^2
f! = get_equations_of_motion(H,2)

function online_eval_energy(θ,par,tol,max_rep,N)
    s = Series(Mean(), Variance())
    for n in 1:max_rep
        OnlineStats.fit!(s,energyMonteCarlo(θ,par,H,2,f!,N,callback=disc_caustic_callback))
        val,σ2 = value.(s.stats)
        if √σ2/val < tol && n>1
            break
        end
    end

    result = value.(s.stats)
    result[1],√result[2]
end

function online_eval_heat(θ,par,tol,max_rep,N)
    s = Series(Mean(), Variance())
    for n in 1:max_rep
        OnlineStats.fit!(s,heatMonteCarlo(θ,par,H,2,f!,N,callback=disc_caustic_callback))
        val,σ2 = value.(s.stats)
        if √σ2/val < tol && n>1
            break
        end
    end

    result = value.(s.stats)
    result[1],√result[2]
end
##

##
μ = .5
θ_min = .1
θ_max = 4
N = 16
scatter_θs = LinRange(θ_min,θ_max,N)
line_θs = LinRange(θ_min,θ_max,4N)

xs = LinRange(-4.5,4.5,160)
ys = LinRange(-4,5,160)
Es,ψs = solveSchrodinger(xs,ys,V;nev=70,par=μ)
Us_ex = [exact_energy(θ,Es) for θ in line_θs]

Us_sc = @showprogress [online_eval_energy(θ,μ,.02,10,3*10^5) for θ in scatter_θs]
Us_c = classical_energy(line_θs,2,μ,potential=V,reltol=1e-3)

#jldsave("Results/Nelson/energy_μ=$μ.jld2";scatter_θs,line_θs,μ,Us_ex,Us_sc,Us_c)
##
using Plots,LaTeXStrings

default(label=false,width=3,size=(354,250), markersize = 4, msw=0, 
palette=:Set1_3, tickfontsize=8,
guidefontsize=8, fontfamily="Computer Modern",dpi=1000,grid=false)

p = plot(line_θs,Us_ex,
        ylabel=L"E/\hbar \omega",xlabel=L"\omega \theta",
        annotations = ((.9,.9), Plots.text(L"\mu=%$μ",10)),
        yaxis=:log) 

        plot!(scatter_θs,first.(Us_sc),label=false,markershape=:diamond,ribbon=last.(Us_sc),linewidth=2)
plot!(line_θs,Us_c,label=false,line=:dot)
##
μ = 2
θ_min = .05
θ_max = 2
N = 16
scatter_θs = LinRange(θ_min,θ_max,N)
line_θs = LinRange(θ_min,θ_max,4N)

xs = LinRange(-4.5,4.5,160)
ys = LinRange(-4,5,160)
Es,ψs = solveSchrodinger(xs,ys,V;nev=70,par=μ)
Cs_ex = [exact_heat(θ,Es) for θ in line_θs]

online_eval_heat(2,μ,1e-2,10,10^5)

Cs_sc = @showprogress [online_eval_heat(θ,μ,.01,60,4*10^5) for θ in scatter_θs]
Cs_c = ThreadsX.map(θ -> classical_heat(θ,2,μ,potential=V,integration_limit=350*exp(-2θ),alg=CubatureJLp(),reltol=1e-9,abstol=1e-8),line_θs)
##
using Plots,LaTeXStrings

default()
default(label=false,width=3,size=(354,250), markersize = 4, msw=0, 
palette=:Set1_3, tickfontsize=8,
guidefontsize=8, fontfamily="Computer Modern",dpi=1000,grid=false,framestyle = :box)

p = plot(line_θs,Cs_ex,
        ylabel=L"c/k",xlabel=L"\omega \theta",
        annotations = ((.9,.8), Plots.text(L"\mu=%$μ",10)),ylims=(0,2.1),xlims=(0,θ_max)) 

plot!(scatter_θs,first.(Cs_sc),label=false,markershape=:diamond,linewidth=2,ribbon=last.(Cs_sc))
plot!(line_θs,Cs_c,label=false,line=:dot)
##
jldsave("Results/Nelson/heat_μ=$μ.jld2";scatter_θs,line_θs,μ,Cs_ex,Cs_sc,Cs_c)
#png("Plots/Nelson/μ=$μ")
##

using Plots,LaTeXStrings

default()
default(label=false,width=3,size=(354,250), markersize = 4, msw=0, 
palette=:Set1_3, tickfontsize=8,
guidefontsize=8, fontfamily="Computer Modern",dpi=1000,grid=false,framestyle = :box)

function make_plot(line_θs,scatter_θs,exact,sc,classical,μ,show_tick_legend)
    p = plot(line_θs,exact,
        ylabel=L"\log_2 E",
        annotations = ((.9,.9), Plots.text(L"\mu=%$μ",10)),
        yaxis=:log,ylims=(0.9*last(classical),1.1*first(classical)),
        yticks=(2 .^ (-1:.8:4 ),-1:.8:4))
    
    if show_tick_legend
        plot!(p,xlabel=L"\theta",bottom_margin=-50Plots.mm)
    else
        plot!(p,xformatter=_->"",bottom_margin=-7.3Plots.mm)
    end    

    plot!(scatter_θs,first.(sc),label=false,
    markershape=:diamond,ribbon=last.(sc),linewidth=0)
    plot!(line_θs,classical,label=false,line=:dot)
end

ps = []

for μ ∈ (.5,1,2)
    _load(name) = load("Results/Nelson/energy_μ=$μ.jld2",name)
    push!(ps,make_plot(_load("line_θs"),_load("scatter_θs"),
    _load("Us_ex"),_load("Us_sc"),_load("Us_c"),μ,length(ps)==2))
end

plot(ps[1],ps[2],ps[3],size=(354,600),layout=(4,1),left_margin=3.5Plots.mm)
#png("Plots/Nelson/energies_nelson")
##
using Plots,LaTeXStrings

default()
default(label=false,width=3,size=(354,250), markersize = 4, msw=0, 
palette=:Set1_3, tickfontsize=8,
guidefontsize=8, fontfamily="Computer Modern",dpi=1000,grid=false,framestyle = :box)

function make_plot(line_θs,scatter_θs,exact,sc,classical,μ,show_tick_legend)
    p = plot(line_θs,exact,
        ylabel=L"c/k",
        annotations = ((.9,.75), Plots.text(L"\mu=%$μ",10)),
        )
    
    if show_tick_legend
        plot!(p,xlabel=L"\theta",bottom_margin=-50Plots.mm)
    else
        plot!(p,xformatter=_->"",bottom_margin=-7.3Plots.mm)
    end    

    plot!(scatter_θs,first.(sc),label=false,
    markershape=:diamond,ribbon=last.(sc),linewidth=0)
    plot!(line_θs,classical,label=false,line=:dot)
end

ps = []

for μ ∈ (.5,1,2)
    _load(name) = load("Results/Nelson/heat_μ=$μ.jld2",name)
    push!(ps,make_plot(_load("line_θs"),_load("scatter_θs"),
    _load("Cs_ex"),_load("Cs_sc"),_load("Cs_c"),μ,length(ps)==2))
end

plot(ps[1],ps[2],ps[3],size=(354,600),layout=(4,1),left_margin=3.5Plots.mm)
#png("Plots/Nelson/heats_nelson")