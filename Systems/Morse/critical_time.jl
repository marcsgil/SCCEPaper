using Plots,LaTeXStrings

default(label=false,width=3,dpi=1000,size=(354,250), markersize = 5, msw=0, 
palette=:Dark2_3, tickfontsize=8, labelfontsize=12,
guidefontsize=8, fontfamily="Computer Modern",grid=false,framestyle = :box)

θc(ϵ) = log(1/√ϵ+√(1/ϵ-1))/√(1-ϵ)

ϵs = LinRange(.001,1,1024)
p = plot(ϵs,θc.(ϵs),xlabel=L"\epsilon",ylabel=L"\omega s_c",
xlims=(0,1),ylims=(.9,3.1),size=(354,200),dpi=300)
##
png(p,"Plots/Morse/critical_time.png")
