using Plots,LaTeXStrings

default(label=false,width=3,dpi=100,size=(354,250), markersize = 5, msw=0, 
palette=:Dark2_3, tickfontsize=10, 
guidefontsize=12, fontfamily="Computer Modern")

θc(ϵ) = log(1/√ϵ+√(1/ϵ-1))/√(1-ϵ)

ϵs = LinRange(.001,1,1024)
p = plot(ϵs,θc.(ϵs),xlabel=L"\epsilon",ylabel=L"\omega \theta_c",
xlims=(0,1),ylims=(.9,3.1),guidefontsize=15,size=(354,200),dpi=300)
##
png(p,"Plots/critical_time2.png")