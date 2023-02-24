using SCCanonicalEnsemble
using ClassicalCanonicalEnsemble
using ExactCanonicalEnsemble
##
H(x,χ) = sum(abs2,x)/2+χ*(sum(abs2,x)/2)^2
Es(χ) = [ (n+1/2) + χ*(n+1/2)^2 for n in 0:10^5 ]

sc_energy(θs::AbstractArray,χ) = energy_NF(θs,Polynomial([0,1,χ]),2000.)
sc_heat(θs::AbstractArray,χ) = heat_NF(θs,Polynomial([0,1,χ]),2000.)
##
θs = LinRange(.2,10,128)
χ = .5
Us_ex =  exact_energy(θs,χ,Es)
Us_sc = sc_energy(θs,χ)

Us_c = classical_energy(θs,H,1,χ)
##
C_ex(θs,χ)
Us_sc = C_sc(θs,χ)
Cs_c = classical_heat.(θs,H,1,χ)
##
using Plots,LaTeXStrings

default(label=false,width=3,size=(354,250), markersize = 5, msw=0, 
palette=:Set1_3, tickfontsize=8,
guidefontsize=8, fontfamily="Computer Modern",dpi=1000)

function make_plot(θs,exact,sc,classical,χ,show_tick_legend)
    p = plot(θs,exact,
        ylabel=L"U/\hbar \omega_0",
        annotations = ((.9,.9), Plots.text(L"\chi=%$χ",10)))
    
    if show_tick_legend
        plot!(p,xlabel=L"\omega_0 \theta",bottom_margin=-3Plots.mm)
    else
        plot!(p,xformatter=_->"",bottom_margin=-7.5Plots.mm)
    end    

    plot!(θs,sc,label=false,line=:dash)
    plot!(θs,classical,label=false,line=:dot)
end

θs = LinRange(.5,10,128)
ps = []

for χ ∈ [.1,.3,.6,1]
    Us_ex = exact_energy(θs,χ,Es)
    Us_sc = sc_energy(θs,χ)
    Us_c = classical_energy(θs,H,1,χ)
    push!(ps,make_plot(θs,Us_ex,Us_sc,Us_c,χ,length(ps)==3))
end

plot(ps[1],ps[2],ps[3],ps[4],size=(354,700),layout=(4,1),left_margin=3Plots.mm)
#png("Plots/Kerr/energies_kerr")
##
using Plots,LaTeXStrings

default(label=false,width=3,size=(354,250), markersize = 5, msw=0, 
palette=:Set1_3, tickfontsize=8,
guidefontsize=8, fontfamily="Computer Modern",dpi=1000)

function make_plot(θs,exact,sc,classical,χ,show_tick_legend)
    p = plot(θs,exact,
        ylabel=L"c/k",
        annotations = ((.9,.5), Plots.text(L"\chi=%$χ",10)))
    
    if show_tick_legend
        plot!(p,xlabel=L"\omega_0 \theta",bottom_margin=-3Plots.mm)
    else
        plot!(p,xformatter=_->"",bottom_margin=-7.8Plots.mm)
    end    

    plot!(θs,sc,label=false,line=:dash)
    plot!(θs,classical,label=false,line=:dot)
end

θs = LinRange(.1,10,128)
ps = []

for χ ∈ [.1,.3,.6,1]
    Cs_ex = exact_heat(θs,χ,Es)
    Cs_sc = sc_heat(θs,χ)
    Cs_c = classical_heat(θs,H,1,χ)
    push!(ps,make_plot(θs,Cs_ex,Cs_sc,Cs_c,χ,length(ps)==3))
end

plot(ps[1],ps[2],ps[3],ps[4],size=(354,700),layout=(4,1),left_margin=3Plots.mm,ylims=(-.2,.99))
#png("Plots/Kerr/heats_kerr")