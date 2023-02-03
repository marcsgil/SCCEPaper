using DifferentialEquations,StaticArrays,ForwardDiff,Symbolics,LinearAlgebra

X = [-.01,-log(2)+.01]
tspan = (0,1.)
χ = .12
##
function f(u,χ,t)
    SA[-4χ*u[3],
    ( -exp(-u[4])*cos(0.5*u[1]) + exp(-2*u[4])*cos(u[1]) )/χ,
    (exp(-u[4])*sin(0.5*u[1]) - exp(-2*u[4])*sin(u[1]) )/(2χ),
    -χ*u[2],
    (exp(-u[4])*sin(0.5*u[1]) - exp(-2*u[4])*sin(u[1]) )/(2χ)*u[1] - χ*u[2]*u[2]]
end

prob = ODEProblem( f, SVector{5}(vcat(zero(X),X,zeros(eltype(X)))), tspan, χ)
solve(prob)
@benchmark solve(prob,save_everystep=false,save_start=false)
##
function get_sol(x,prob)
    _prob = remake(prob,u0=SVector{5}(vcat(zero(x),x,zeros(eltype(x)))))
    sols = solve(_prob,save_everystep=false,save_start=false)
    view(sols,3:4,length(sols))
end

get_sol(X,prob)
ForwardDiff.jacobian(x->get_sol(x,prob),X)

@benchmark ForwardDiff.jacobian(x->get_sol(x,prob),X)
##



myexp(x) = exp(real(x))*cis(imag(x))
H(x,χ) = χ*x[1]^2+(1-myexp(-x[2]))^2/(4χ)

g = get_equations_of_motion(H,1)

u0 = @SVector rand(5)
g(u0,χ,0)
f(u0,χ,0)



prob2 = ODEProblem( g, SVector{5}(vcat(zero(X),X,zeros(eltype(X)))), tspan, χ)
solve(prob2)
@benchmark solve(prob2,save_everystep=false,save_start=false)

get_sol(X,prob2)
ForwardDiff.jacobian(x->get_sol(x,prob2),X)
@benchmark ForwardDiff.jacobian(x->get_sol(x,prob2),X)