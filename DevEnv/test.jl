using DifferentialEquations,SparseDiffTools,LinearAlgebra
using Integrals, IntegralsCuba,IntegralsCubature

V(q,χ) = χ*q^4/4+(1/2-χ)*q^2
H(x,χ) = x[1]^2/2 + V(x[2],χ)
#f! = get_equations_of_motion(H,1,"χ")


function f!(du,u,χ)
    du[1] = -2*u[3]
    du[2] = -4u[4]*(1/2 - χ) - χ*(4u[4]*(u[4]^2 - u[1]^2/4) - 2*u[4]*u[1]^2)/2
    du[3] = χ*(-2*u[1]*u[4]^2 - u[1]*(u[4]^2 - u[1]^2/4))/2 - u[1]*(1/2 - χ)
    du[4] = -u[2]/2
    nothing
end

function u0(X)
    #Given the initial center X, builds the initial conditions for the initial value problem.
    T = eltype(X)
    N = length(X)
    vcat(zero(X),X,vec(vcat(zeros(eltype(X),N,N),I(N))),zero(T))
end

function phase_space_dim(u)
    Int((-1+√(2length(u)-1))/2)
end

function extract_jac_x(u)
    N = phase_space_dim(u)
    view(reshape((@view u[2N+1:end-1]),2N,N),N+1:2N,:)
end

extract_det_jac(u) = det(extract_jac_x(u))

function F!(du,u,par,f!,J,N)
    #Differential equation for y and x
    f!(view(du,1:2N),view(u,1:2N),par)

    #Differential equation for the jacobians
    J.op.u .= view(u,1:2N)

    for j in 1:N
        mul!(view(du, 2j*N+1:2*(j+1)*N ),J,view(u,2j*N+1:2*(j+1)*N ))
    end

    #Differential equation for the area Δ
    du[end] = view(u,1:N) ⋅ view(du,N+1:2N)

    nothing
end

function annul!(integrator)
    integrator.u = zero(integrator.u)
end

disc_caustic_cross_contidion(u,t,integrator) = extract_det_jac(u) < 0 && extract_det_jac(integrator.uprev) > 0
disc_caustic_callback = DiscreteCallback(disc_caustic_cross_contidion,annul!,save_positions=(false,false))

function energy_integrand(Xs,(θ,par))
    u₀ = u0(view(Xs,:,1))
    N = phase_space_dim(u₀)

    J = JacVec((du,u) -> f!(du,u,par),view(u₀,1:2N))

    prob = ODEProblem((du,u,par,t)->F!(du,u,par,J,N),u₀,(0,θ/2),par)

    function prob_func(prob,i,repeat)
        remake(prob,u0=u0(view(Xs,:,i)))
    end

    function output_func(sol,i)
        Z = exp( sol[end,1] - θ*H(view(Xs,:,i),par) )*√abs(extract_det_jac(sol.u[1])) 
        output = [Z,Z*H(sol[N+1:2N,1],par)]
        isfinite(sum(output)) ? (output,false) : (zero(output),false)
    end

    ensemble_prob = EnsembleProblem(prob;prob_func,output_func)

    solve(ensemble_prob,BS3(),trajectories=size(Xs,2),reltol=1e-2,abstol=1e-3,
    save_everystep=false,save_start=false,callback=disc_caustic_callback,verbose=false) |> stack
end
##
ub = fill(5.,2)
lb = -ub
χ = 1.
θ = 3.
energy_integrand(rand(2,16),(θ,χ))

prob = IntegralProblem(energy_integrand, lb, ub,(θ,χ),nout = 2, batch = 64)
U = solve(prob, CubatureJLh(), reltol = 1e-2, abstol = 1e-3)
@benchmark solve(prob, CubatureJLh(), reltol = 1e-3, abstol = 1e-3)

U[2] / U[1]
##
using SolveSchrodinger

xs = LinRange(-10,10,2048)
Es,ψs = solveSchrodinger(xs,V;par=χ)
ex_U(θ,Es) = sum(E->E*exp(-θ*E),Es)/sum(E->exp(-θ*E),Es)

ex_U(θ,Es)