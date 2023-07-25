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
    solve(prob, Tsit5(), reltol = 1e-6, abstol = 1e-6, save_start=false, save_everystep=false)[end-length(X)+1:end,1]
end

X = rand(2)
flow(X,0,1.)
@benchmark flow($X,1.,1.)

χ = 1.
t = 1.
x = rand(2)
prob = NonlinearProblem((X,par) -> flow(X,t,par) - x,rand(2),χ)
X = solve(prob)
@benchmark solve(prob)

flow(X,t,χ)
##

