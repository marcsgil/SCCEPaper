using SCCanonicalEnsemble, DiffEqPhysics, NonlinearSolve, FastGaussQuadrature, SpecialFunctions,ClassicalOrthogonalPolynomials, WignerTransform, JLD2, UnPack
includet("../../plot_config.jl")
ξ(χ,q) = exp(-q)/χ
s(χ,n) = 1/χ-2n-1
normalization(χ,n) = s(χ,n)*gamma(n+1)/gamma(1/χ-n)
N_max(χ) = Int(floor(.5*(1/χ-1)))
E(n,χ) = (n+.5)-χ*(n+.5)^2

function ψ(y,χ,n) 
	result = √(  normalization(χ,n)*ξ(χ,y)^s(χ,n)*exp(-ξ(χ,y))  )*laguerrel(n,s(χ,n),ξ(χ,y))
	ifelse(isnan(result),zero(result),result)
end

H(x,χ) = χ*x[1]^2+(1-exp(-x[2]))^2/(4χ)

function f!(du,u,χ)
    du[1] = -4χ*u[3]
    du[2] = ( -exp(-u[4])*cos(0.5*u[1]) + exp(-2*u[4])*cos(u[1]) )/χ
    du[3] = (exp(-u[4])*sin(0.5*u[1]) - exp(-2*u[4])*sin(u[1]) )/(2χ)
    du[4] = -χ*u[2]
end

function dy!(Δy,y,x,χ,t)
    Δy[1] = -4χ*x[1]
    Δy[2] = ( -exp(-x[2])*cos(0.5*y[1]) + exp(-2*x[2])*cos(y[1]) )/χ
end

function dx!(Δx,y,x,χ,t)
    Δx[1] = (exp(-x[2])*sin(0.5*y[1]) - exp(-2*x[2])*sin(y[1]) )/(2χ)
    Δx[2] = -χ*y[2]
end

sol_output(sol, i, Xs, θs, par, H, H2) = sol[size(Xs,1)+1:2*size(Xs,1),1],false

function flow(X,θ,par)
    integrand(reshape(X,length(X),1), [θ], par, f!, H, sol_output; save_start=false, save_everystep=false, verbose=false) |> vec
end

function getNodesAndWeights(χ,N=300)
    Ps,wPs = gausslegendre(N)
    Qs,wQs = gausschebyshev(N,3)

    [[√(1-Q^2)*P/(2χ),-log(1-Q)] for P in Ps, Q in Qs] |> vec |> stack, [w1*w2 for w1 in wPs, w2 in wQs] |> vec
end
##
χ = .01
θ = 3.
ps = LinRange(-15,15,512)
qs = LinRange(-0.3,.3,512)
xs = [[p,q] for q ∈ qs, p ∈ ps]
Xs = similar(xs)

Threads.@threads for n ∈ eachindex(xs)
    prob = NonlinearProblem((X,par) -> flow(X,θ,par) - xs[n],xs[n],χ)
    Xs[n] = solve(prob, reltol = 1e-9)
end

Ws = integrand(Xs |> vec |> stack, [θ], χ, f!, H, partition_output_quadrature; save_start=false, save_everystep=false, verbose=false, reltol=1e-8,abstol=1e-8)
Ws = reshape(Ws,length(qs),length(ps)) / (sum(Ws) * (qs[2] - qs[1]) * (ps[2] - ps[1]))
Ps = dropdims(sum(Ws,dims=1),dims=1) * (qs[2] - qs[1])
Qs = dropdims(sum(Ws,dims=2),dims=2) * (ps[2] - ps[1])

quantum_Ws = sum(n -> wigner_transform(q -> ψ(q,χ,n),qs,ps) * exp(-θ*E(n-1,χ)) |> permutedims , 0:N_max(χ)) / sum(n -> exp(-θ*E(n-1,χ)), 0:N_max(χ))
quantum_Ps = dropdims(sum(quantum_Ws,dims=1),dims=1) * (qs[2] - qs[1])
quantum_Qs = dropdims(sum(quantum_Ws,dims=2),dims=2) * (ps[2] - ps[1])

classical_Z(θ,χ) = dawson(√(θ/4χ))/√θ
classical_Ws = [exp(-θ * H([p,q],χ)) for q ∈ qs, p ∈ ps]
classical_Ws /= (sum(classical_Ws) * (qs[2] - qs[1]) * (ps[2] - ps[1]))
sum(classical_Ws) * (qs[2] - qs[1]) * (ps[2] - ps[1])

classical_Ps = dropdims(sum(classical_Ws,dims=1),dims=1) * (qs[2] - qs[1])
classical_Qs = dropdims(sum(classical_Ws,dims=2),dims=2) * (ps[2] - ps[1])
##
file =  jldopen("/home/marcsgil/Code/SCCanonicalEnsemble/Results/Morse/wigner.jld2")
@unpack ps,qs,Ps,Qs,Ws,quantum_Ps,quantum_Qs,quantum_Ws,classical_Ps,classical_Qs,classical_Ws = file
##
fig = Figure(resolution=(1000, 600))
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[2, 1])
ax3 = Axis(fig[2, 2])
 
lines!(ax1, qs , quantum_Qs; color=colors[1], label = "Quantum")
lines!(ax3, quantum_Ps, ps; color=colors[1])

lines!(ax1, qs , Qs; color=colors[2], linestyle = :dot, linewidth=8, label = "SC")
heatmap!(ax2, qs, ps, Ws, colorlimits = (0, maximum(Ws)), colormap = :plasma)
lines!(ax3, Ps, ps; color=colors[2], linestyle = :dot, linewidth=8)

lines!(ax1, qs , classical_Qs; color=colors[3], linestyle=:dash, label="Classical")
lines!(ax3, classical_Ps, ps; color=colors[3],  linestyle=:dash)

Colorbar(fig[:, end+1], limits = (0, maximum(Ws)), colormap = :plasma, label = L"W(p,q)")

hideydecorations!(ax3,)
hidexdecorations!(ax1,)

leg = Legend(fig[1, 2], ax1)
leg.tellheight = true

colsize!(fig.layout, 1, Relative(3 / 4))
rowsize!(fig.layout, 1, Relative(1 / 4))
colgap!(fig.layout, 10)
rowgap!(fig.layout, 10)

ax2.xlabel = L"q"
ax2.ylabel = L"p"

ax1.ylabel = L"W(q)"
ax3.xlabel = L"W(p)"



fig
#save("/home/marcsgil/Code/SCCanonicalEnsemble/Plots/Morse/wigner.svg",fig)
##
#jldsave("/home/marcsgil/Code/SCCanonicalEnsemble/Results/Morse/wigner.jld2";ps,qs,Ps,Qs,Ws,quantum_Ps,quantum_Qs,quantum_Ws,classical_Ps,classical_Qs,classical_Ws)