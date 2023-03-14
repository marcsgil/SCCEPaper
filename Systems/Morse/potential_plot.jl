using Plots,LaTeXStrings

default(label=false,width=3,dpi=1000,size=(354,250), markersize = 5, msw=0, 
palette=:Dark2_3, tickfontsize=8, labelfontsize=10,
guidefontsize=8, fontfamily="Computer Modern",grid=false,framestyle = :box)

V(q) = (1-exp(-q))^2

qs = LinRange(-.8,4,1024)
p = plot(qs,V.(qs),xlabel=L"a(r-r_e)",ylabel=L"V",size=(354,200),dpi=300,
yticks=([0,.5,1,1.5],[L"0",L"D/2",L"D",L"3D/2"]))
##
png(p,"Plots/Morse/morse_potential.png")