using SCCanonicalEnsemble
using ClassicalCanonicalEnsemble
##
Es(χ) = [ (n+1/2) + χ*(n+1/2)^2 for n in 0:10^5 ]

function U_ex(θ::Number,Es::AbstractArray)
    sum(E->E*exp(-θ*E),Es)/sum(E->exp(-θ*E),Es)
end

function U_ex(θs::AbstractArray,χ::Number)
    es = Es(χ)
    map(θ->U_ex(θ,es),θs)
end

function C_ex(θ::Number,Es::AbstractArray)
    Z = sum(E->exp(-θ*E),Es)
    U = sum(E->E*exp(-θ*E),Es)/Z
    U2 = sum(E->E^2*exp(-θ*E),Es)/Z
    θ^2*( U2 - U^2 )
end

function C_ex(θs::AbstractArray,χ::Number)
    es = Es(χ)
    map(θ->C_ex(θ,es),θs)
end


U_sc(θs::AbstractArray,χ) = energy_NF(θs,Polynomial([0,1,χ]),2000.)
C_sc(θs::AbstractArray,χ) = heat_NF(θs,Polynomial([0,1,χ]),2000.)

function U_sc(θs::AbstractArray,χs::AbstractArray)
    reduce(hcat,map(χ->U(θs,χ),χs))
end
##
θs = LinRange(.2,10,128)
χ = .5
Us_ex =  U_ex(θs,χ)
Us_sc = U_sc(θs,χ)
H(x,χ) = sum(abs2,x)/2+χ*(sum(abs2,x)/2)^2
Us_c = classical_energy.(θs,H,1,χ)
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
    Us_ex = U_ex(θs,χ)
    Us_sc = U_sc(θs,χ)
    Us_c = classical_energy.(θs,H,1,χ)
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
    Cs_ex = C_ex(θs,χ)
    Cs_sc = C_sc(θs,χ)
    Cs_c = classical_heat.(θs,H,1,χ)
    push!(ps,make_plot(θs,Cs_ex,Cs_sc,Cs_c,χ,length(ps)==3))
end

plot(ps[1],ps[2],ps[3],ps[4],size=(354,700),layout=(4,1),left_margin=3Plots.mm,ylims=(-.2,.99))
png("Plots/Kerr/heats_kerr")