using Plots,LaTeXStrings

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
bottom_margin = -5Plots.mm,
left_margin = -5Plots.mm,
dpi=1000
)
png("Plots/Nelson/nelson_potential")