using FFTW,Interpolations

function get_wigner_grid(qs,χs)
	[ (q+0.5*χ,q-0.5*χ) for q in qs, χ in χs ]
end

function wigner_from_func(ψ,qs,ps::DFTGrid)
	wigner_grid = get_wigner_grid(qs,reciprocal_grid(ps))
	auto_correlation = map( Q -> ψ(Q[1])*conj(ψ(Q[2])),wigner_grid )

	real(fftshift(fft(ifftshift(auto_correlation,2),2),2))/(2*ps.max_val)
end

function wigner_from_vec(ψs::AbstractArray,qs,ps::DFTGrid)
    wigner_from_func(cubic_spline_interpolation(qs, ψs,extrapolation_bc=0),ps,qs)
end

function calculate_many_wigners(ψs::AbstractArray,qs,ps::DFTGrid,max::Int)
	#ψs should be a matrix whose columns are the wavefunctions
	wigners = zeros(ComplexF64,(length(qs),ps.N,min(size(ψs,2),max)))
	wigner_grid = get_wigner_grid(qs,reciprocal_grid(ps))
	plan = plan_fft(wigners[:,:,1],1)
	for n in axes(wigners,3)
		ψ = cubic_spline_interpolation(qs, view(ψs,:, n),extrapolation_bc=0)
		auto_correlation = map( Q -> ψ(Q[1])*conj(ψ(Q[2])),wigner_grid )
		wigners[:,:,n] = fftshift(fft(ifftshift(auto_correlation,2),2),2)
	end
	real(wigners)/(2*ps.max_val)
end

function calculate_many_wigners(ψ::Function,qs,ps::DFTGrid,max::Int)
	wigners = zeros(ComplexF64,(length(qs),ps.N,max))
	wigner_grid = get_wigner_grid(qs,reciprocal_grid(ps))
	plan = plan_fft(wigners[:,:,1],1)
	for n in axes(wigners,3)
		auto_correlation = map( Q -> ψ(Q[1],n-1)*conj(ψ(Q[2],n-1)),wigner_grid )
		wigners[:,:,n] = fftshift(fft(ifftshift(auto_correlation,2),2),2)
	end
	real(wigners)/(2*ps.max_val)
end