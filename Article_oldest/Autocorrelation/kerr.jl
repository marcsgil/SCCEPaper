using FFTW
using Integrals,IntegralsCubature
using CairoMakie

ψ₀(q) = π^(-.25)*exp(-.5*q^2)

function ψ_cs(q::Real,α)
    Q = √2*real(α)
    P = √2*imag(α)

    ψ₀(q-Q)*cis(P*(q-.5Q))
end

function ψ_kerr(q,a::Int,b::Int,α)
    N = mod(b,4) == 0 ? b÷2 : b
    
    c = map(n->cis(-2π*a*(n-1)^2/b ),1:N) |> ifft

    sum(k-> c[k]*ψ_cs.(q, α*cis(-2π*( (k-1)/N + a/b ) ) ), 1:N)
end

function auto_correlation(a::Int,b::Int,α)
    prob = IntegralProblem((q,α)->conj(ψ_kerr(q,0,1,α))ψ_kerr(q,a,b,α),-Inf,Inf,α)
    solve(prob,HCubatureJL()).u
end

function W_cs(x,Q,P)
    exp(-( (x[1]-Q)^2+(x[2]-P)^2 ))/π
end

function sc_auto_correlation_integrand(X,(t,Q,P))
    r2 = X[1]^2+X[2]^2
    ϕ = t
    s,c = sincos(ϕ)
    det_jac = c^2-2ϕ*s*c
    cis(( t*r2-s )*r2*.5)*W_cs(c*X,Q,P)√abs(det_jac)
    W_cs(cos(t)*X,Q,P)
end

function sc_auto_correlation(t,α)
    prob = IntegralProblem(sc_auto_correlation_integrand,fill(-Inf,2),fill(Inf,2),(t,α))
    solve(prob,HCubatureJL())
end
##
using BenchmarkTools
α = 5/√2
t = .01
@benchmark W_cs(1,1,$α)
@benchmark sc_auto_correlation_integrand([1,1],(.01,1,1))
@code_warntype sc_auto_correlation_integrand([1,1],(t,α))
##
grid = LinRange(-50,50,1024)
zs = [sc_auto_correlation_integrand([X[1],X[2]],(1.5,5,0)) |> abs2 for X in Iterators.product(grid,grid)]
heatmap(grid,grid,zs,colormap=cgrad([:red,:white,:blue]),colorrange=(-1/2π,1/2π))

##
α = 5/√2
b = 4096
as = 0:512
ts = as.*(π/(2b))
cs = auto_correlation.(as,b,α)
sc_auto_correlation(.05,α)
##
fig,ax=lines(ts,abs2.(cs))
ylims!(ax,-.01,.3)
fig