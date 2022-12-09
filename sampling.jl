using CairoMakie
using Roots,Distributions

function sample(f,N,support,y_max)
    result = Vector{Float64}(undef,N)

    Δ = support[2] - support[1]
    
    n = 1
    while n ≤ N
        x = support[1] + Δ*rand()
        if f(x) ≥ y_max*rand()
            result[n] = x
            n += 1
        end
    end
    
    result
end

function sample(f,N,support)
   max = map(f,LinRange(minimum(support),maximum(support),2048)) |> maximum
   sample(f,N,support,1.1*max)
end

function sample(f,N)
    xmax = find_zero(x->f(x)-1e-4,1) |> abs
    sample(f,N,(-xmax,xmax))
end

function validate_samples(f,samples)
    xmax = find_zero(x->f(x)-1e-4,1)

    fig,ax = hist(samples,normalization = :pdf,bins=100)
    xs = LinRange(-xmax,xmax,1024)

    lines!(ax,xs,f.(xs)/(sum(f.(xs))*abs(xs[2]-xs[1])),color=:red)
    fig
end
##
#=V(q,χ) = χ*q^4/4+(1/2-χ)*q^2
θ = 6
χ = 1
f(x) = exp(-θ*V(x,χ))

samples = sample(f,10^6)
validate_samples(f,samples)=#