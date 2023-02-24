module ExactCanonicalEnsemble

export exact_energy,exact_heat

function exact_energy(θ::Real,Es::AbstractArray)
    sum(E->E*exp(-θ*E),Es)/sum(E->exp(-θ*E),Es)
end

function exact_energy(θs::AbstractArray,Es::AbstractArray)
    map(θ->exact_energy(θ,Es),θs)
end

function exact_energy(θ,par,Es::Function)
    exact_energy(θ,Es(par))
end

function exact_heat(θ::Real,Es::AbstractArray)
    Z = sum(E->exp(-θ*E),Es)
    U = sum(E->E*exp(-θ*E),Es)/Z
    U2 = sum(E->E^2*exp(-θ*E),Es)/Z
    θ^2*( U2 - U^2 )
end

function exact_heat(θs::AbstractArray,Es::AbstractArray)
    map(θ->exact_heat(θ,Es),θs)
end

function exact_heat(θ,par,Es::Function)
    exact_heat(θ,Es(par))
end

end
