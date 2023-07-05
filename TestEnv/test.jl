using AdvancedMH,Distributions,LinearAlgebra,MCMCChains,StatsPlots

# Use a struct instead of `typeof(density)` for sake of readability.

density(x) = -(x[1]^2/2 + x[2]^4 - x[2]^2)

m1 = DensityModel(density)
p1 = MetropolisHastings(RandomWalkProposal(MvNormal(zeros(2),I)))
c1 = sample(m1, p1, 10^5; param_names = [:p,:q], chain_type=Chains)

propertynames(c1)
propertynames(c1.value)
c1.value.data
size(c1)

histogram(c1.value.data[:,1,1],normalize=true)
##
@benchmark sample(m1, p1, 10^5; param_names = [:p,:q], chain_type=Chains)

function metropolisHastings(f,N)
  samples = Vector{Float64}(undef,N)

  samples[1] = randn()

  for n in 2:length(samples)
      x₀ = samples[n-1]
      x = randn() + x₀

      if f(x)/f(x₀) ≥ rand()
          samples[n] = x
      else
          samples[n] = x₀
      end
  end
  
  samples
end

function metropolisHastings(f,N,d)
  samples = Array{Float64}(undef,d,N)

  samples[:,1] = randn(d)

  for n in 2:size(samples,2)
      x₀ = view(samples,:,n-1)
      x = randn(d) .+ x₀

      if f(x)/f(x₀) ≥ rand()
          samples[:,n] = x
      else
          samples[:,n] = x₀
      end
  end
  
  samples
end


pdf = exp ∘ density
c2 = metropolisHastings(pdf,10^5,2)
histogram(c2,normalize=true)
@benchmark metropolisHastings($pdf,10^5,2)