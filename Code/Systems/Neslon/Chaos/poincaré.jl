using DynamicalSystems, StaticArrays
using OrdinaryDiffEq: Vern9

using CairoMakie, MathTeXEngine, ColorSchemes

colors = colorschemes[:hsv]
theme = Theme(
    fontsize=32,
    fonts = (; regular = texfont()),
    Axis=(xlabelsize=36, xlabelpadding=0,
    ylabelsize=36, ylabelpadding=0,
    xticklabelsize = 28, yticklabelsize=28,
    xgridvisible=false, ygridvisible=false,
    xtickalign=1, ytickalign=1,
    yticksize=12, xticksize=12),
    Lines = (linewidth = 5,),
    Legend = (labelsize=28,)
)

set_theme!()
set_theme!(theme)

diffeq = (alg=Vern9(), abstol=1e-9, reltol=1e-9)

function nelson(u, μ, t)
    du1 = u[3]
    du2 = u[4]
    du3 = 2u[1] * (u[2] - u[1]^2/2 - μ)
    du4 = u[1]^2 - 2u[2]
    return SVector(du1, du2, du3, du4)
end

u = rand(4)
nelson(u,μ,0) ≈ nelson2(u, μ, 0)


##
E = 4.8
μ = 2.0
u0s = Vector{Float64}[]

for _ ∈ 1:length(colors)
    for _ ∈ 1:1000
        global x = randn()
        global y = randn()
        global K = E - V(x, y, μ)
        if K > 0
            break
        end
    end
    px = (rand() - 0.5) * 2√(2K)
    py = √(2K - px^2)
    push!(u0s,[x,y,px,py])
end
##
fig = Figure()
ax = Axis(fig[1,1])
ax.xlabel = L"x"
ax.ylabel = L"p_x"

for (n,u0) ∈ enumerate(u0s)
    ds = CoupledODEs(nelson, u0, μ; diffeq)
    data = poincaresos(ds, (2,0.),1e4)
    xs = data[:,1]
    pxs = data[:,3]
    scatter!(ax,xs,pxs,markersize=3,color = colors[n],alpha=0.2)
end

fig
##
save("/home/marcsgil/Code/SCCanonicalEnsemble/Plots/Nelson/poincare.pdf",fig)