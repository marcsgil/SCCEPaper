using Elliptic.Jacobi,GLMakie,Roots
import Elliptic.Jacobi.sn,Elliptic.Jacobi.cn,Elliptic.Jacobi.dn

sn(u) = sn(u,1/2)
cn(u) = cn(u,1/2)
dn(u) = dn(u,1/2)

f(u) = sn(u)/dn(u)
##
E(x) = (x[2]^2+x[1]^4)/2
A(x) = (2E(x))^(1/4)
ω₀(x) = (8E(x))^(1/4)
ϕ(x) = find_zero(ϕ->f(ϕ)-√2*x[1]/A(x),0) 

function x(t::Number,x₀)
    freq = ω₀(x₀)
    u = freq*t+ϕ(x₀)
    amp = A(x₀)/√2
    func = f(u)
    
    (amp*func,amp*freq*cn(u)*(1+func^2/2))
end

function x(ts::AbstractArray,x₀)
    result = Array{eltype(x₀)}(undef,2,length(ts))
    freq = ω₀(x₀)
    phase = ϕ(x₀)
    amp = A(x₀)/√2

    for n in axes(result,2)
        u = freq*ts[n]+phase
        func = f(u)
        result[1,n] = amp*func
        result[2,n] = amp*freq*cn(u)*(1+func^2/2)
    end
    
    result
end
##
ts = LinRange(0,6,100)
x₀ = rand(2)

xs = x(ts,x₀)

fig,ax = lines(ts,@view xs[1,:])
lines!(ax,ts,@view xs[2,:])
fig
#[ps[1],qs[1]] ≈ x₀


#lines(ts,[E([ps[n],qs[n]]) for n in eachindex(ts)])
#map(E,eachcol(xs))