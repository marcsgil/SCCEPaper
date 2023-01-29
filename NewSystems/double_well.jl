include("../sc_solution.jl")
include("../sc_functions.jl")
include("../ex_functions.jl")
include("../plotting_functions.jl")

V(q,χ) = χ*q^4/4+(1/2-χ)*q^2
H(x,χ) = x[1]^2/2 + V(x[2],χ)
fy,fx = get_equations_of_motion(H,1)
##
χ = 1.
xs = LinRange(-10,10,2048)
Es,ψs = solveSchrodinger(xs,V;par=χ)

θ = 3.
ex_U(θ,Es)



energyMonteCarlo(θ,χ,H,1,fy,fx,10^5)
##
N = 64
ps = LinRange(-12,12,N)
qs = LinRange(-6,6,N)
xs = [SA[p,q] for p in ps, q in qs]

θs = LinRange(0,3,64)
χ=2

sols = solve_equations(θs,χ,fy,fx,(par)->(xs,ones(N^2)),H,
abstol=1e-12,reltol=1e-12,alg=Vern8(),output_func=analysis_output,callback=custom_callback)
##
using GLMakie
GLMakie.activate!()
##
fig = Figure(resolution=(1600,700), fontsize=24)
ax1 = GLMakie.Axis(fig[1,1],xlabel="q",ylabel="p",title="Integrando F. de Partição",xticks=-5:1:5)
ax2 = GLMakie.Axis(fig[1,3],xlabel="q",ylabel="p",title="Determinante Jacobiano",xticks=-5:1:5)

sl = Slider(fig[2, :], range = eachindex(θs), startvalue = 1)

title = Label(fig[0, :], L"\theta = %$(round(θs[1],digits=2)),\chi=%$χ", fontsize = 30)
lift(sl.value) do n
    title.text = L"\theta = %$(round(θs[n],digits=2)),\chi=%$χ"
end

Zs = lift(n->reshape(view(sols,1,n,:),N,N) |> transpose,sl.value)
dets = lift(n->reshape(view(sols,2,n,:),N,N) |> transpose,sl.value)

hm1 = heatmap!(ax1,qs,ps,Zs,colormap=:afmhot,colorrange=(0,3))
Colorbar(fig[1, 2], hm1)

hm2 = heatmap!(ax2,qs,ps,dets,colormap=:redsblues,colorrange=(-1,1))
Colorbar(fig[1, 4], hm2)

fig
##
χ = 1
θ_min = .5
θ_max = 3
N = 16
θs_sc = LinRange(θ_min,θ_max,N)
θs_ex = LinRange(θ_min,θ_max,4N)

xs = LinRange(-10,10,2048)
Es,ψs = solveSchrodinger(xs,V;par=χ)

Us_ex = [ex_U(θ,Es) for θ in θs_ex]
Us_sc = @showprogress [energyMonteCarlo(θ,χ,H,1,fy,fx,10^6) for θ in θs_sc]

comparison_plot(θs_ex,θs_sc,Us_ex,Us_sc,L"\chi = %$χ")