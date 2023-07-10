using Plots,LaTeXStrings

default()

V(x,y,μ) = (x^2/2-y)^2 + μ*x^2

V(x,y) = V(x,y,2)

x = LinRange(-4,4,512)
y = LinRange(-4,13,512)
surface(x, y, V,
camera = (75 , 40),
xlabel=L"x",
ylabel=L"y",
zlabel=L"V(x,y)",
color=:gnuplot2,
size=(400,420),colorbar=false,
top_margin = -6Plots.mm,
bottom_margin = -10Plots.mm,
left_margin = -12Plots.mm,
dpi=1000,
guidefontsize=16,
yticks = -3:3:13,
tickfontsize=10
)
#png("Plots/Nelson/nelson_potential")