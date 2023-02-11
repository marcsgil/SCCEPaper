using SCCanonicalEnsemble
using ClassicalCanonicalEnsemble
##
Es(χ) = [ (n+1/2) + χ*(n+1/2)^2 for n in 0:10^4 ]
U_ex(θ,χ) = sum(E->E*exp(-θ*E),Es(χ))/sum(E->exp(-θ*E),Es(χ))

U_sc(θs::AbstractArray,χ) = energy_NF(θs,Polynomial([0,1,χ]),2000.)

function U_sc(θs::AbstractArray,χs::AbstractArray)
    reduce(hcat,map(χ->U(θs,χ),χs))
end
##
θs = LinRange(.2,10,128)
χ = .1

Us_ex = U_ex.(θs,χ)
Us_sc = U_sc(θs,χ)
H(x,χ) = sum(abs2,x)/2+χ*(sum(abs2,x)/2)^2
Us_c = classical_energy.(θs,H,1,χ)
##
using Plots,LaTeXStrings

default(label=false,width=3,size=(354,250), markersize = 5, msw=0, 
palette=:Set1_3, tickfontsize=8,
guidefontsize=8, fontfamily="Computer Modern")

function make_plot(Us_ex,Us_sc,Us_c,θs,χ,show_tick_legend)
    if show_tick_legend
        plot(θs,Us_ex,ylabel=L"U",xlabel=L"\theta",
        annotations = ((.9,.9), Plots.text(L"\chi=%$χ",10)),bottom_margin=-3Plots.mm)
    else
        plot(θs,Us_ex,xformatter=_->"",ylabel=L"U",
        annotations = ((.9,.9), Plots.text(L"\chi=%$χ",10)),bottom_margin=-7.5Plots.mm)
    end

    plot!(θs,Us_sc,label=false,line=:dash)
    plot!(θs,Us_c,label=false,line=:dot)
end

θs = LinRange(.5,10,128)
ps = []

for χ ∈ [.1,.3,.6,1]
    Us_ex = U_ex.(θs,χ)
    Us_sc = U_sc(θs,χ)
    Us_c = classical_energy.(θs,H,1,χ)
    push!(ps,make_plot(Us_ex,Us_sc,Us_c,θs,χ,length(ps)==3))
end

#plot!(ps[4],xlabel=L"\theta")

plot(ps[1],ps[2],ps[3],ps[4],size=(354,700),layout=(4,1),left_margin=3Plots.mm)
#png("Plots/Kerr/energies_kerr")