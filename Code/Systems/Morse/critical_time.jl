using Plots,LaTeXStrings

default()
default(label=false,width=4, markersize = 4, msw=0, 
palette=:Set1_3, tickfontsize=14,
guidefontsize=16, fontfamily="Computer Modern",dpi=1000,grid=false,framestyle = :box)

θc(ϵ) = log(1/√ϵ+√(1/ϵ-1))/√(1-ϵ)

ϵs = LinRange(.001,1,1024)
p = plot(ϵs,θc.(ϵs),xlabel=L"\epsilon",ylabel=L"\omega s_c",
xlims=(0,1),ylims=(.9,3.1),size=(600,400),dpi=300,left_margin=2Plots.mm,bottom_margin=4Plots.mm)
##
png(p,"Plots/Morse/critical_time.png")
