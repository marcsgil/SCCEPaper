module WaveToWigner

using FFTW,Interpolations
using CUDA,CUDA.CUFFT

include("dft_utils.jl")

export waveToWinger,calculate_many_wigners

function waveToWinger(ψ_interp,ps,qs)
    #ψ_interp = cubic_spline_interpolation(qs, ψ, extrapolation_bc=0)
	auto_correlation = map( Q -> ψ_interp(Q[1]+Q[2]/2)*conj(ψ_interp(Q[1]-Q[2]/2)), Iterators.product(qs,reciprocal_grid(ps)) )

	real(fftshift(fft(ifftshift(auto_correlation,2),2),2))/(2*last(ps))
end

function waveToWinger(ψs::AbstractArray,ps,qs)
    waveToWinger(cubic_spline_interpolation(qs, ψs, extrapolation_bc=0),ps,qs)
end

function calculate_many_wigners(ψs,ps,qs,max=Inf)
	#ψs should be a matrix whose columns are the wavefunctions
	wigners = Array{complex(eltype(ψs))}(undef, size(qs)..., size(ps)..., min(size(ψs)[end],max))
	grid = Iterators.product(qs...,reciprocal_grid(ps)...)
	plan = plan_fft( first(eachslice(wigners,dims=ndims(wigners))) )
	for n in axes( wigners, ndims(wigners) )
		auto_correlation = map( (q,χ) -> ψ(q+χ/2)*conj(ψ(q-χ/2)), Iterators.product(qs...,reciprocal_grid(ps)...) )
		wigners[ntuple(i->:,ndims(wigners)-1),n] = fftshift(fft(ifftshift(auto_correlation,2),2),2)
	end
	real(wigners)/(2*last(ps))
end

end