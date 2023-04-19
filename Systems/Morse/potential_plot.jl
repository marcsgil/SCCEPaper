using Plots,LaTeXStrings

default()
default(label=false,width=4, markersize = 4, msw=0, 
palette=:Set1_3, tickfontsize=14,
guidefontsize=16, fontfamily="Computer Modern",dpi=1000,grid=false,framestyle = :box)

V(q) = (1-exp(-q))^2

qs = LinRange(-.8,4,1024)
p = plot(qs,V.(qs),xlabel=L"a(r-r_e)",ylabel=L"V",dpi=300,
yticks=([0,.5,1,1.5],[L"0",L"D/2",L"D",L"3D/2"]))
##
png(p,"Plots/Morse/morse_potential.png")