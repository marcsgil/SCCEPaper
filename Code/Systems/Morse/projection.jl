ps_1 = (2,-2,-8,8)
ps_2 = (2,-22/3,8,-8/3)

function φ(ξ,ps_1,ps_2)
    if abs(ξ) > 1
        zero(ξ)
    elseif abs(ξ) < 0.5
        evalpoly(abs(ξ),ps_1)
    else
        evalpoly(abs(ξ),ps_2)
    end
end

function φ(ξ,ps_1,ps_2,ϵ)
    prod(ξ-> φ(ξ/ϵ,ps_1,ps_2), ξ) / ϵ^(ndims(ξ))
end


ξs = LinRange(-2,2,512)
@benchmark φ(1.,$ps_1,$ps_2)

ξs = rand(10)
@benchmark φ($ξs,$ps_1,$ps_2,1)

lines(ξs,[φ(ξ,ps_1,ps_2) for ξ ∈ ξs])
##

using SCCanonicalEnsemble, NonlinearSolve
includet("../../plot_config.jl")
##
H(x,χ) = χ*x[1]^2+(1-exp(-x[2]))^2/(4χ)

function f!(du,u,χ,t)
    du[1] = -4χ*u[3]
    du[2] = ( -exp(-u[4])*cos(0.5*u[1]) + exp(-2*u[4])*cos(u[1]) )/χ
    du[3] = (exp(-u[4])*sin(0.5*u[1]) - exp(-2*u[4])*sin(u[1]) )/(2χ)
    du[4] = -χ*u[2]
end

function flow(X,t,par)
    prob = ODEProblem{true}(f!,vcat(zero(X),X),(0,t),par)
    solve(prob, BS3(), reltol = 1e-2, abstol = 1e-3, save_start=false, save_everystep=false)[end-length(X)+1:end,1]
end
##
χ = .12
θ = 1.5
ps = LinRange(-5,5,128)
qs = LinRange(-1,2,128)
grid = [[p,q] for q ∈ qs, p ∈ ps]
Xs = similar(grid)

Threads.@threads for n ∈ eachindex(grid)
    prob = NonlinearProblem((X,par) -> flow(X,θ/2,par) - grid[n],grid[n],χ)
    Xs[n] = solve(prob)
end

Ws = integrand(Xs |> vec |> stack, θ, χ, f!, H, wigner_output; save_start=false, save_everystep=false, verbose=false)
heatmap(qs,ps,reshape(Ws,128,128),colormap=:hot)